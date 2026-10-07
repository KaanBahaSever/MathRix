#pragma once

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <ostream>
#include <vector>

namespace mathrix {

/// Polynomial with real coefficients, stored in ascending order:
/// p(x) = c[0] + c[1] x + c[2] x^2 + ...
class Polynomial {
public:
    Polynomial() = default;
    Polynomial(std::initializer_list<double> ascending) : c_(ascending) { trim(); }
    explicit Polynomial(std::vector<double> ascending) : c_(std::move(ascending)) { trim(); }

    /// Degree (-1 for the zero polynomial).
    int degree() const { return static_cast<int>(c_.size()) - 1; }
    const std::vector<double>& coefficients() const { return c_; }
    double coefficient(size_t power) const { return power < c_.size() ? c_[power] : 0.0; }

    /// Horner evaluation.
    double operator()(double x) const {
        double r = 0.0;
        for (size_t i = c_.size(); i-- > 0;) r = r * x + c_[i];
        return r;
    }

    Polynomial derivative() const {
        if (c_.size() <= 1) return {};
        std::vector<double> d(c_.size() - 1);
        for (size_t i = 1; i < c_.size(); ++i) d[i - 1] = c_[i] * static_cast<double>(i);
        return Polynomial(d);
    }

    /// Antiderivative with the given integration constant.
    Polynomial integral(double constant = 0.0) const {
        std::vector<double> d(c_.size() + 1, 0.0);
        d[0] = constant;
        for (size_t i = 0; i < c_.size(); ++i) d[i + 1] = c_[i] / static_cast<double>(i + 1);
        return Polynomial(d);
    }

    /// Definite integral over [a, b].
    double integrate(double a, double b) const {
        const Polynomial P = integral();
        return P(b) - P(a);
    }

    Polynomial operator+(const Polynomial& o) const {
        std::vector<double> r(std::max(c_.size(), o.c_.size()), 0.0);
        for (size_t i = 0; i < c_.size(); ++i) r[i] += c_[i];
        for (size_t i = 0; i < o.c_.size(); ++i) r[i] += o.c_[i];
        return Polynomial(r);
    }
    Polynomial operator-(const Polynomial& o) const { return *this + o * -1.0; }
    Polynomial operator*(double s) const {
        std::vector<double> r = c_;
        for (double& x : r) x *= s;
        return Polynomial(r);
    }
    friend Polynomial operator*(double s, const Polynomial& p) { return p * s; }
    Polynomial operator*(const Polynomial& o) const {
        if (c_.empty() || o.c_.empty()) return {};
        std::vector<double> r(c_.size() + o.c_.size() - 1, 0.0);
        for (size_t i = 0; i < c_.size(); ++i)
            for (size_t j = 0; j < o.c_.size(); ++j) r[i + j] += c_[i] * o.c_[j];
        return Polynomial(r);
    }
    bool operator==(const Polynomial& o) const { return c_ == o.c_; }

    /// Real root in [a, b] by bisection + Newton (requires a sign change).
    double rootIn(double a, double b, double tol = 1e-12) const {
        double fa = (*this)(a);
        if (fa == 0.0) return a;
        const Polynomial d = derivative();
        for (int it = 0; it < 200 && b - a > tol; ++it) {
            double m = 0.5 * (a + b);
            const double dm = d(m);
            if (dm != 0.0) {
                const double n = m - (*this)(m) / dm;
                if (n > a && n < b) m = n;
            }
            const double fm = (*this)(m);
            if (fm == 0.0) return m;
            if ((fm > 0) == (fa > 0)) {
                a = m;
                fa = fm;
            } else {
                b = m;
            }
        }
        return 0.5 * (a + b);
    }

private:
    std::vector<double> c_;
    void trim() {
        while (!c_.empty() && c_.back() == 0.0) c_.pop_back();
    }
};

inline std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
    if (p.degree() < 0) return os << "0";
    bool first = true;
    for (int i = p.degree(); i >= 0; --i) {
        const double c = p.coefficient(static_cast<size_t>(i));
        if (c == 0.0) continue;
        os << (first ? (c < 0 ? "-" : "") : (c < 0 ? " - " : " + "));
        const double a = c < 0 ? -c : c;
        if (a != 1.0 || i == 0) os << a;
        if (i >= 1) os << "x";
        if (i >= 2) os << "^" << i;
        first = false;
    }
    return os;
}

}  // namespace mathrix
