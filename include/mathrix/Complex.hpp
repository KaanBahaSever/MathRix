#pragma once

#include <cmath>
#include <complex>
#include <ostream>
#include <vector>

namespace mathrix {

/// Complex number z = re + i im (double precision).
/// Interoperates with std::complex<double> in both directions.
struct Complex {
    double re = 0.0;
    double im = 0.0;

    constexpr Complex() = default;
    constexpr Complex(double real, double imag = 0.0) : re(real), im(imag) {}
    Complex(const std::complex<double>& c) : re(c.real()), im(c.imag()) {}

    /// r (cos t + i sin t)
    static Complex polar(double r, double theta) { return {r * std::cos(theta), r * std::sin(theta)}; }
    static constexpr Complex i() { return {0.0, 1.0}; }

    constexpr double real() const { return re; }
    constexpr double imag() const { return im; }
    std::complex<double> toStd() const { return {re, im}; }

    constexpr Complex conj() const { return {re, -im}; }
    /// |z|
    double abs() const { return std::hypot(re, im); }
    /// |z|^2
    constexpr double norm() const { return re * re + im * im; }
    /// arg z in (-pi, pi]
    double arg() const { return std::atan2(im, re); }

    constexpr Complex operator-() const { return {-re, -im}; }
    constexpr Complex operator+(const Complex& o) const { return {re + o.re, im + o.im}; }
    constexpr Complex operator-(const Complex& o) const { return {re - o.re, im - o.im}; }
    constexpr Complex operator*(const Complex& o) const { return {re * o.re - im * o.im, re * o.im + im * o.re}; }
    /// Division with Smith's algorithm (avoids overflow for large components).
    Complex operator/(const Complex& o) const {
        if (std::abs(o.re) >= std::abs(o.im)) {
            const double r = o.im / o.re, d = o.re + o.im * r;
            return {(re + im * r) / d, (im - re * r) / d};
        }
        const double r = o.re / o.im, d = o.re * r + o.im;
        return {(re * r + im) / d, (im * r - re) / d};
    }
    Complex& operator+=(const Complex& o) { return *this = *this + o; }
    Complex& operator-=(const Complex& o) { return *this = *this - o; }
    Complex& operator*=(const Complex& o) { return *this = *this * o; }
    Complex& operator/=(const Complex& o) { return *this = *this / o; }

    constexpr bool operator==(const Complex& o) const { return re == o.re && im == o.im; }
    constexpr bool operator!=(const Complex& o) const { return !(*this == o); }
};

// Mixed real/complex arithmetic (the Complex(double) constructor covers `z op x`).
constexpr Complex operator+(double a, const Complex& z) { return Complex(a) + z; }
constexpr Complex operator-(double a, const Complex& z) { return Complex(a) - z; }
constexpr Complex operator*(double a, const Complex& z) { return Complex(a) * z; }
inline Complex operator/(double a, const Complex& z) { return Complex(a) / z; }

inline double abs(const Complex& z) { return z.abs(); }
inline double arg(const Complex& z) { return z.arg(); }
constexpr Complex conj(const Complex& z) { return z.conj(); }

/// e^z
inline Complex exp(const Complex& z) { return Complex::polar(std::exp(z.re), z.im); }
/// Principal natural logarithm.
inline Complex log(const Complex& z) { return {std::log(z.abs()), z.arg()}; }
/// Principal square root.
inline Complex sqrt(const Complex& z) { return Complex(std::sqrt(z.toStd())); }
/// Principal power z^w.
inline Complex pow(const Complex& z, const Complex& w) {
    if (z.re == 0.0 && z.im == 0.0) return (w.re == 0.0 && w.im == 0.0) ? Complex(1.0) : Complex(0.0);
    return exp(w * log(z));
}
inline Complex pow(const Complex& z, double x) { return pow(z, Complex(x)); }
/// Integer power by repeated squaring (exact for small n).
inline Complex pow(Complex z, int n) {
    if (n < 0) return 1.0 / pow(z, -n);
    Complex r(1.0);
    while (n) {
        if (n & 1) r *= z;
        z *= z;
        n >>= 1;
    }
    return r;
}
inline Complex sin(const Complex& z) { return {std::sin(z.re) * std::cosh(z.im), std::cos(z.re) * std::sinh(z.im)}; }
inline Complex cos(const Complex& z) { return {std::cos(z.re) * std::cosh(z.im), -std::sin(z.re) * std::sinh(z.im)}; }

/// All n distinct n-th roots of z, w_k = |z|^(1/n) e^{i (arg z + 2 pi k) / n}, k = 0..n-1.
inline std::vector<Complex> nthRoots(const Complex& z, int n) {
    std::vector<Complex> out;
    if (n <= 0) return out;
    const double r = std::pow(z.abs(), 1.0 / n);
    const double pi = 3.14159265358979323846;
    for (int k = 0; k < n; ++k) out.push_back(Complex::polar(r, (z.arg() + 2.0 * pi * k) / n));
    return out;
}

/// Sets components smaller than `eps` (relative to |z|) to exactly zero, removing round-off
/// residue such as e^{i pi} = -1 + 1.2e-16 i.
inline Complex chop(const Complex& z, double eps = 1e-12) {
    const double scale = eps * (1.0 + z.abs());
    return {std::abs(z.re) < scale ? 0.0 : z.re, std::abs(z.im) < scale ? 0.0 : z.im};
}

/// True when |a - b| <= tol.
inline bool approxEqual(const Complex& a, const Complex& b, double tol = 1e-12) { return (a - b).abs() <= tol; }

inline std::ostream& operator<<(std::ostream& os, const Complex& z) {
    if (z.im == 0.0) return os << z.re;
    if (z.re == 0.0) return os << z.im << "i";
    return os << z.re << (z.im < 0 ? " - " : " + ") << std::abs(z.im) << "i";
}

namespace literals {
/// 3.0_i  ->  Complex(0, 3)
constexpr Complex operator""_i(long double x) { return {0.0, static_cast<double>(x)}; }
constexpr Complex operator""_i(unsigned long long x) { return {0.0, static_cast<double>(x)}; }
}  // namespace literals

}  // namespace mathrix
