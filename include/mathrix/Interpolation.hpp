#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mathrix {

/// Interpolation schemes for 1-D tabulated data.
///  - Step:        zero-order hold (value of the left sample)
///  - Linear:      piecewise linear
///  - CubicSpline: natural C2 cubic spline (smooth, may overshoot)
///  - Akima:       Akima spline (smooth, robust against outliers)
///  - Pchip:       monotone piecewise cubic Hermite (no overshoot; best for thrust curves)
enum class Interp { Step, Linear, CubicSpline, Akima, Pchip };

/// What happens outside of the sampled range.
enum class Extrapolation { Clamp, Zero, Linear };

inline const char* toString(Interp m) {
    switch (m) {
        case Interp::Step: return "step";
        case Interp::Linear: return "linear";
        case Interp::CubicSpline: return "spline";
        case Interp::Akima: return "akima";
        case Interp::Pchip: return "pchip";
    }
    return "linear";
}

inline Interp interpFromString(const std::string& s) {
    if (s == "step" || s == "zoh") return Interp::Step;
    if (s == "linear" || s == "lineer") return Interp::Linear;
    if (s == "spline" || s == "cubic" || s == "cubicspline" || s == "natural") return Interp::CubicSpline;
    if (s == "akima") return Interp::Akima;
    if (s == "pchip" || s == "monotone" || s == "hermite") return Interp::Pchip;
    throw std::invalid_argument("unknown interpolation method '" + s + "'");
}

inline const char* toString(Extrapolation e) {
    switch (e) {
        case Extrapolation::Clamp: return "clamp";
        case Extrapolation::Zero: return "zero";
        case Extrapolation::Linear: return "linear";
    }
    return "clamp";
}

inline Extrapolation extrapolationFromString(const std::string& s) {
    if (s == "clamp" || s == "hold") return Extrapolation::Clamp;
    if (s == "zero") return Extrapolation::Zero;
    if (s == "linear") return Extrapolation::Linear;
    throw std::invalid_argument("unknown extrapolation mode '" + s + "'");
}

/// 1-D interpolator. All cubic schemes are stored in Hermite form (values + slopes), so
/// evaluation, derivative and exact integration share one code path.
class Interpolator1D {
public:
    Interpolator1D() = default;

    Interpolator1D(std::vector<double> x, std::vector<double> y,
                   Interp method = Interp::Linear,
                   Extrapolation extrapolation = Extrapolation::Clamp)
        : x_(std::move(x)), y_(std::move(y)), method_(method), extrap_(extrapolation) {
        if (x_.size() != y_.size()) throw std::invalid_argument("Interpolator1D: x and y sizes differ");
        if (x_.empty()) throw std::invalid_argument("Interpolator1D: no samples");
        sortAndMerge();
        computeSlopes();
    }

    bool empty() const { return x_.empty(); }
    size_t size() const { return x_.size(); }
    const std::vector<double>& xs() const { return x_; }
    const std::vector<double>& ys() const { return y_; }
    Interp method() const { return method_; }
    Extrapolation extrapolation() const { return extrap_; }
    double xMin() const { return x_.front(); }
    double xMax() const { return x_.back(); }

    double operator()(double x) const { return eval(x); }

    double eval(double x) const {
        if (x_.empty()) return 0.0;
        if (x_.size() == 1) return extrap_ == Extrapolation::Zero && x != x_[0] ? 0.0 : y_[0];
        if (x < x_.front()) return extrapolate(x, 0);
        if (x > x_.back()) return extrapolate(x, x_.size() - 1);
        const size_t i = segment(x);
        return evalSegment(i, x);
    }

    double derivative(double x) const {
        if (x_.size() < 2) return 0.0;
        if (x < x_.front() || x > x_.back()) {
            if (extrap_ != Extrapolation::Linear) return 0.0;
            return x < x_.front() ? d_.front() : d_.back();
        }
        const size_t i = segment(x);
        const double h = x_[i + 1] - x_[i];
        if (method_ == Interp::Step) return 0.0;
        if (method_ == Interp::Linear) return (y_[i + 1] - y_[i]) / h;
        const double t = (x - x_[i]) / h;
        const double dh00 = 6 * t * t - 6 * t, dh10 = 3 * t * t - 4 * t + 1;
        const double dh01 = -6 * t * t + 6 * t, dh11 = 3 * t * t - 2 * t;
        return (dh00 * y_[i] + dh01 * y_[i + 1]) / h + dh10 * d_[i] + dh11 * d_[i + 1];
    }

    /// Exact integral of the interpolant over [a, b] (inside the sampled range; outside
    /// regions follow the extrapolation rule).
    double integrate(double a, double b) const {
        if (x_.size() < 2 || a == b) return 0.0;
        if (a > b) return -integrate(b, a);
        double total = 0.0;
        // Left extrapolated region.
        if (a < x_.front()) {
            const double e = std::min(b, x_.front());
            total += simpson(a, e);
            a = e;
            if (a >= b) return total;
        }
        // Right extrapolated region.
        double tail = 0.0;
        if (b > x_.back()) {
            const double s = std::max(a, x_.back());
            tail = simpson(s, b);
            b = s;
            if (a >= b) return total + tail;
        }
        size_t i = segment(a);
        while (a < b && i + 1 < x_.size()) {
            const double e = std::min(b, x_[i + 1]);
            total += segmentIntegral(i, a, e);
            a = e;
            ++i;
        }
        return total + tail;
    }

private:
    std::vector<double> x_, y_, d_;
    Interp method_ = Interp::Linear;
    Extrapolation extrap_ = Extrapolation::Clamp;

    void sortAndMerge() {
        std::vector<std::pair<double, double>> p(x_.size());
        for (size_t i = 0; i < x_.size(); ++i) p[i] = {x_[i], y_[i]};
        std::stable_sort(p.begin(), p.end(), [](auto& a, auto& b) { return a.first < b.first; });
        x_.clear();
        y_.clear();
        for (auto& [px, py] : p) {
            if (!x_.empty() && std::abs(px - x_.back()) < 1e-12 * (1.0 + std::abs(px))) {
                y_.back() = py;  // duplicate abscissa: keep last value
                continue;
            }
            x_.push_back(px);
            y_.push_back(py);
        }
    }

    size_t segment(double x) const {
        auto it = std::upper_bound(x_.begin(), x_.end(), x);
        size_t i = static_cast<size_t>(std::distance(x_.begin(), it));
        if (i == 0) return 0;
        i -= 1;
        if (i >= x_.size() - 1) i = x_.size() - 2;
        return i;
    }

    double evalSegment(size_t i, double x) const {
        const double h = x_[i + 1] - x_[i];
        const double t = (x - x_[i]) / h;
        switch (method_) {
            case Interp::Step: return (x >= x_[i + 1]) ? y_[i + 1] : y_[i];
            case Interp::Linear: return y_[i] + (y_[i + 1] - y_[i]) * t;
            default: break;
        }
        const double t2 = t * t, t3 = t2 * t;
        const double h00 = 2 * t3 - 3 * t2 + 1, h10 = t3 - 2 * t2 + t;
        const double h01 = -2 * t3 + 3 * t2, h11 = t3 - t2;
        return h00 * y_[i] + h10 * h * d_[i] + h01 * y_[i + 1] + h11 * h * d_[i + 1];
    }

    double extrapolate(double x, size_t end) const {
        switch (extrap_) {
            case Extrapolation::Clamp: return y_[end];
            case Extrapolation::Zero: return 0.0;
            case Extrapolation::Linear: {
                double slope;
                if (method_ == Interp::Step) slope = 0.0;
                else if (method_ == Interp::Linear) {
                    slope = end == 0 ? (y_[1] - y_[0]) / (x_[1] - x_[0])
                                     : (y_[end] - y_[end - 1]) / (x_[end] - x_[end - 1]);
                } else slope = d_[end];
                return y_[end] + slope * (x - x_[end]);
            }
        }
        return y_[end];
    }

    double segmentIntegral(size_t i, double a, double b) const {
        if (method_ == Interp::Step) return y_[i] * (b - a);
        if (method_ == Interp::Linear) return 0.5 * (evalSegment(i, a) + evalSegment(i, b)) * (b - a);
        const double m = 0.5 * (a + b);  // Simpson is exact for cubics
        return (b - a) / 6.0 * (evalSegment(i, a) + 4.0 * evalSegment(i, m) + evalSegment(i, b));
    }

    double simpson(double a, double b) const {
        const double m = 0.5 * (a + b);
        return (b - a) / 6.0 * (eval(a) + 4.0 * eval(m) + eval(b));
    }

    void computeSlopes() {
        const size_t n = x_.size();
        d_.assign(n, 0.0);
        if (n < 2) return;
        std::vector<double> h(n - 1), m(n - 1);
        for (size_t i = 0; i + 1 < n; ++i) {
            h[i] = x_[i + 1] - x_[i];
            m[i] = (y_[i + 1] - y_[i]) / h[i];
        }
        if (n == 2 || method_ == Interp::Linear || method_ == Interp::Step) {
            for (size_t i = 0; i < n; ++i) d_[i] = m[std::min(i, n - 2)];
            return;
        }
        switch (method_) {
            case Interp::CubicSpline: naturalSpline(h); break;
            case Interp::Akima: akima(m); break;
            case Interp::Pchip: pchip(h, m); break;
            default: break;
        }
    }

    void naturalSpline(const std::vector<double>& h) {
        const size_t n = x_.size();
        std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), r(n, 0.0);
        b[0] = 2.0 / h[0];
        c[0] = 1.0 / h[0];
        r[0] = 3.0 * (y_[1] - y_[0]) / (h[0] * h[0]);
        for (size_t i = 1; i + 1 < n; ++i) {
            a[i] = 1.0 / h[i - 1];
            b[i] = 2.0 * (1.0 / h[i - 1] + 1.0 / h[i]);
            c[i] = 1.0 / h[i];
            r[i] = 3.0 * ((y_[i] - y_[i - 1]) / (h[i - 1] * h[i - 1]) + (y_[i + 1] - y_[i]) / (h[i] * h[i]));
        }
        a[n - 1] = 1.0 / h[n - 2];
        b[n - 1] = 2.0 / h[n - 2];
        r[n - 1] = 3.0 * (y_[n - 1] - y_[n - 2]) / (h[n - 2] * h[n - 2]);
        // Thomas algorithm.
        for (size_t i = 1; i < n; ++i) {
            const double w = a[i] / b[i - 1];
            b[i] -= w * c[i - 1];
            r[i] -= w * r[i - 1];
        }
        d_[n - 1] = r[n - 1] / b[n - 1];
        for (size_t i = n - 1; i-- > 0;) d_[i] = (r[i] - c[i] * d_[i + 1]) / b[i];
    }

    void akima(const std::vector<double>& m) {
        const size_t n = x_.size();
        if (n < 3) return;
        // Extended slope array: index k in ext corresponds to m[k-2].
        std::vector<double> e(n + 3);
        for (size_t i = 0; i < n - 1; ++i) e[i + 2] = m[i];
        e[1] = 2.0 * e[2] - e[3];
        e[0] = 2.0 * e[1] - e[2];
        e[n + 1] = 2.0 * e[n] - e[n - 1];
        e[n + 2] = 2.0 * e[n + 1] - e[n];
        for (size_t i = 0; i < n; ++i) {
            const double m0 = e[i], m1 = e[i + 1], m2 = e[i + 2], m3 = e[i + 3];
            const double w1 = std::abs(m3 - m2), w2 = std::abs(m1 - m0);
            d_[i] = (w1 + w2 < 1e-12) ? 0.5 * (m1 + m2) : (w1 * m1 + w2 * m2) / (w1 + w2);
        }
    }

    void pchip(const std::vector<double>& h, const std::vector<double>& m) {
        const size_t n = x_.size();
        for (size_t k = 1; k + 1 < n; ++k) {
            if (m[k - 1] * m[k] <= 0.0) {
                d_[k] = 0.0;
            } else {
                const double w1 = 2.0 * h[k] + h[k - 1];
                const double w2 = h[k] + 2.0 * h[k - 1];
                d_[k] = (w1 + w2) / (w1 / m[k - 1] + w2 / m[k]);
            }
        }
        auto endSlope = [](double h0, double h1, double m0, double m1) {
            double d = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
            if ((d > 0) != (m0 > 0) || m0 == 0.0) d = 0.0;
            else if (((m0 > 0) != (m1 > 0)) && std::abs(d) > 3.0 * std::abs(m0)) d = 3.0 * m0;
            return d;
        };
        d_[0] = endSlope(h[0], h[1], m[0], m[1]);
        d_[n - 1] = endSlope(h[n - 2], h[n - 3], m[n - 2], m[n - 3]);
    }
};

/// 2-D table z(x, y) on a rectangular grid. Interpolates along x with the chosen scheme
/// and linearly along y. Outside the grid values are clamped.
class Table2D {
public:
    Table2D() = default;
    /// `z[j][i]` is the value at (x[i], y[j]).
    Table2D(std::vector<double> x, std::vector<double> y, std::vector<std::vector<double>> z,
            Interp methodX = Interp::Linear)
        : x_(std::move(x)), y_(std::move(y)), z_(std::move(z)), method_(methodX) {
        if (y_.empty() || x_.empty()) throw std::invalid_argument("Table2D: empty grid");
        if (z_.size() != y_.size()) throw std::invalid_argument("Table2D: row count must equal y size");
        for (auto& row : z_) {
            if (row.size() != x_.size()) throw std::invalid_argument("Table2D: column count must equal x size");
            rows_.emplace_back(x_, row, method_, Extrapolation::Clamp);
        }
        if (!std::is_sorted(y_.begin(), y_.end())) throw std::invalid_argument("Table2D: y must be ascending");
    }

    bool empty() const { return rows_.empty(); }
    const std::vector<double>& xs() const { return x_; }
    const std::vector<double>& ys() const { return y_; }
    const std::vector<std::vector<double>>& values() const { return z_; }
    Interp method() const { return method_; }

    double operator()(double x, double y) const {
        if (rows_.empty()) return 0.0;
        if (rows_.size() == 1 || y <= y_.front()) return rows_.front()(x);
        if (y >= y_.back()) return rows_.back()(x);
        auto it = std::upper_bound(y_.begin(), y_.end(), y);
        const size_t j = static_cast<size_t>(std::distance(y_.begin(), it)) - 1;
        const double t = (y - y_[j]) / (y_[j + 1] - y_[j]);
        return rows_[j](x) * (1.0 - t) + rows_[j + 1](x) * t;
    }

private:
    std::vector<double> x_, y_;
    std::vector<std::vector<double>> z_;
    Interp method_ = Interp::Linear;
    std::vector<Interpolator1D> rows_;
};

}  // namespace mathrix
