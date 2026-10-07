#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <ostream>
#include <vector>

#include "Complex.hpp"

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

    /// Horner evaluation at a complex point.
    Complex operator()(const Complex& z) const {
        Complex r;
        for (size_t i = c_.size(); i-- > 0;) r = r * z + c_[i];
        return r;
    }

    /// All roots (real and complex, with multiplicity) by the Durand-Kerner method followed by
    /// Newton polishing. Roots whose imaginary part is negligible are returned as exact reals.
    std::vector<Complex> roots(double tol = 1e-13, int maxIterations = 1000) const {
        const int n = degree();
        std::vector<Complex> z;
        if (n < 1) return z;
        // Monic coefficients.
        std::vector<double> a(c_.size());
        for (size_t i = 0; i < c_.size(); ++i) a[i] = c_[i] / c_.back();
        const Polynomial monic(a);
        // Initial guesses on a circle of the Cauchy root bound.
        double bound = 0.0;
        for (int i = 0; i < n; ++i) bound = std::max(bound, std::abs(a[static_cast<size_t>(i)]));
        bound += 1.0;
        for (int k = 0; k < n; ++k) z.push_back(Complex::polar(bound, 0.4 + 6.283185307179586 * k / n));
        for (int it = 0; it < maxIterations; ++it) {
            double change = 0.0;
            for (int k = 0; k < n; ++k) {
                Complex denom(1.0);
                for (int j = 0; j < n; ++j)
                    if (j != k) denom *= z[static_cast<size_t>(k)] - z[static_cast<size_t>(j)];
                if (denom.norm() == 0.0) denom = Complex(tol, tol);
                const Complex step = monic(z[static_cast<size_t>(k)]) / denom;
                z[static_cast<size_t>(k)] -= step;
                change = std::max(change, step.abs() / (1.0 + z[static_cast<size_t>(k)].abs()));
            }
            if (change < tol) break;
        }
        const Polynomial d = monic.derivative();
        for (auto& r : z) {
            for (int it = 0; it < 3; ++it) {  // Newton polishing (skipped near multiple roots)
                const Complex dv = d(r);
                if (dv.abs() < 1e-8) break;
                r -= monic(r) / dv;
            }
            if (std::abs(r.im) <= 1e-9 * (1.0 + std::abs(r.re))) r.im = 0.0;
            if (std::abs(r.re) <= 1e-12 * (1.0 + std::abs(r.im))) r.re = 0.0;
        }
        std::sort(z.begin(), z.end(), [](const Complex& x, const Complex& y) {
            return x.re != y.re ? x.re < y.re : x.im < y.im;
        });
        return z;
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
