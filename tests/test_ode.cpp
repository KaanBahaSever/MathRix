#include <mathrix/mathrix.hpp>

#include "TestFramework.hpp"

using namespace mathrix;

TEST_CASE("ODE integrators converge with their order") {
    using V = VecN<2>;
    // Harmonic oscillator x'' = -x, exact x = cos t.
    auto f = [](double, const V& y) {
        V d;
        d[0] = y[1];
        d[1] = -y[0];
        return d;
    };
    auto run = [&](Integrator m, double h) {
        V y;
        y[0] = 1.0;
        y[1] = 0.0;
        double t = 0.0;
        const int n = static_cast<int>(std::lround(2.0 / h));
        for (int i = 0; i < n; ++i, t += h) y = fixedStep(m, f, t, y, h);
        return std::abs(y[0] - std::cos(2.0));
    };
    CHECK(run(Integrator::Euler, 0.01) < 0.02);
    CHECK(run(Integrator::Heun, 0.01) < 1e-4);
    CHECK(run(Integrator::RK4, 0.01) < 1e-9);
    const double ratio = run(Integrator::RK4, 0.02) / run(Integrator::RK4, 0.01);
    CHECK(ratio > 12.0 && ratio < 20.0);  // ~2^4
    V y;
    y[0] = 1.0;
    double t = 0.0, h = 0.1;
    while (t < 2.0 - 1e-12) {
        h = std::min(h, 2.0 - t);
        const auto s = dormandPrinceStep(f, t, y, h, 1e-10, 1e-10);
        if (s.error <= 1.0) {
            y = s.y;
            t += h;
        }
        h = nextStepSize(h, s.error);
    }
    CHECK_NEAR(y[0], std::cos(2.0), 1e-8);
}
