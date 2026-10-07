#include <mathrix/mathrix.hpp>

#include "TestFramework.hpp"

using namespace mathrix;
using namespace mathrix::literals;

TEST_CASE("complex: arithmetic and conversions") {
    const Complex a{3, 4}, b{1, -2};
    CHECK(a + b == Complex(4, 2));
    CHECK(a - b == Complex(2, 6));
    CHECK(a * b == Complex(11, -2));
    CHECK(approxEqual(a / b, Complex(-1, 2)));
    CHECK(approxEqual((a / b) * b, a));
    CHECK_NEAR(a.abs(), 5.0, 1e-15);
    CHECK_NEAR(a.norm(), 25.0, 1e-15);
    CHECK(a.conj() == Complex(3, -4));
    CHECK(2.0 * a == Complex(6, 8));
    CHECK(a + 1.0 == Complex(4, 4));
    CHECK(approxEqual(1.0 / Complex(0, 1), Complex(0, -1)));
    CHECK(3.0_i == Complex(0, 3));
    CHECK(1.0 + 2.0_i == Complex(1, 2));
    const std::complex<double> s = a.toStd();
    CHECK(Complex(s * s) == a * a);
    // Smith's division stays accurate for huge components.
    const Complex big{1e300, 1e300};
    CHECK(approxEqual(big / big, Complex(1, 0)));
}

TEST_CASE("complex: polar form and elementary functions") {
    const double pi = kPi;
    CHECK(approxEqual(exp(Complex(0, pi)), Complex(-1, 0)));  // Euler's identity
    CHECK(chop(exp(Complex(0, pi))) == Complex(-1, 0));
    CHECK(approxEqual(Complex::polar(2.0, pi / 2), Complex(0, 2)));
    CHECK_NEAR(Complex(-1, 0).arg(), pi, 1e-15);
    CHECK(approxEqual(sqrt(Complex(-4, 0)), Complex(0, 2)));
    CHECK(approxEqual(log(exp(Complex(0.5, 1.0))), Complex(0.5, 1.0)));
    CHECK(approxEqual(pow(Complex(0, 1), 2), Complex(-1, 0)));
    CHECK(approxEqual(pow(Complex(1, 1), -2), Complex(0, -0.5)));
    CHECK(approxEqual(pow(Complex(0, 1), Complex(0, 1)), Complex(std::exp(-pi / 2), 0)));  // i^i
    const Complex z{0.3, -0.7};
    CHECK(approxEqual(sin(z) * sin(z) + cos(z) * cos(z), Complex(1, 0)));
    const auto r = nthRoots(Complex(1, 0), 4);
    CHECK(r.size() == 4);
    CHECK(approxEqual(r[1], Complex(0, 1)));
    for (const auto& w : r) CHECK(approxEqual(pow(w, 4), Complex(1, 0)));
}

TEST_CASE("polynomial: real and complex roots") {
    const auto r1 = Polynomial{1, 0, 1}.roots();  // x^2 + 1
    CHECK(r1.size() == 2);
    CHECK(approxEqual(r1[0], Complex(0, -1), 1e-12) && approxEqual(r1[1], Complex(0, 1), 1e-12));
    const Polynomial cubic = Polynomial{-1, 1} * Polynomial{-2, 1} * Polynomial{-3, 1};  // (x-1)(x-2)(x-3)
    const auto r2 = cubic.roots();
    for (int k = 0; k < 3; ++k) {
        CHECK_NEAR(r2[static_cast<size_t>(k)].re, k + 1.0, 1e-12);
        CHECK(r2[static_cast<size_t>(k)].im == 0.0);
    }
    const auto r3 = Polynomial{-1, 0, 0, 1}.roots();  // x^3 - 1: cube roots of unity
    CHECK(r3.size() == 3);
    for (const auto& z : r3) CHECK(approxEqual(pow(z, 3), Complex(1, 0), 1e-12));
    const Polynomial q{5, -2, 3, 0, 1, 7};  // degree 5, random coefficients
    for (const auto& z : q.roots()) CHECK(q(z).abs() < 1e-9);
    const auto r4 = Polynomial{1, -2, 1}.roots();  // (x-1)^2, double root
    CHECK_NEAR(r4[0].re, 1.0, 1e-6);
    CHECK_NEAR(r4[1].re, 1.0, 1e-6);
    CHECK(Polynomial{4}.roots().empty());
}
