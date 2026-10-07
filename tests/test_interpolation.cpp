#include <mathrix/mathrix.hpp>

#include "TestFramework.hpp"

using namespace mathrix;

TEST_CASE("interpolation schemes") {
    const std::vector<double> x = {0, 1, 2, 3, 4}, y = {0, 1, 4, 9, 16};
    const Interpolator1D lin(x, y, Interp::Linear);
    CHECK_NEAR(lin(1.5), 2.5, 1e-12);
    CHECK_NEAR(lin(-1), 0.0, 1e-12);  // clamped
    CHECK_NEAR(lin.integrate(0, 4), 22.0, 1e-12);  // trapezoid
    const Interpolator1D step(x, y, Interp::Step);
    CHECK_NEAR(step(1.9), 1.0, 1e-12);
    for (Interp m : {Interp::CubicSpline, Interp::Akima, Interp::Pchip}) {
        const Interpolator1D it(x, y, m);
        for (size_t i = 0; i < x.size(); ++i) CHECK_NEAR(it(x[i]), y[i], 1e-12);
        CHECK_NEAR(it(2.5), 6.25, 0.2);
        // Exact integration of the interpolant (Simpson on cubic segments).
        double num = 0.0;
        for (int k = 0; k < 40000; ++k) num += it((k + 0.5) * 1e-4) * 1e-4;
        CHECK_NEAR(it.integrate(0, 4), num, 1e-6);
    }
    // PCHIP never overshoots monotone data.
    const Interpolator1D p({0, 1, 2, 3}, {0, 0, 1, 1}, Interp::Pchip);
    for (double t = 0; t <= 3; t += 0.01) CHECK(p(t) >= -1e-12 && p(t) <= 1.0 + 1e-12);
    // Unsorted input with duplicates is accepted.
    const Interpolator1D u({2, 0, 1, 1}, {4, 0, 5, 1}, Interp::Linear);
    CHECK_NEAR(u(0.5), 0.5, 1e-12);
    CHECK(interpFromString("spline") == Interp::CubicSpline);
}

TEST_CASE("2-D table") {
    const Table2D t({0.0, 1.0}, {0.0, 1000.0}, {{1.0, 2.0}, {3.0, 4.0}});
    CHECK_NEAR(t(0.5, 0.0), 1.5, 1e-12);
    CHECK_NEAR(t(0.5, 500.0), 2.5, 1e-12);
    CHECK_NEAR(t(2.0, 5000.0), 4.0, 1e-12);
}

