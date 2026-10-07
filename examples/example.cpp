// MathRix in a few lines: rotate a vector, fit a smooth curve, integrate an ODE, solve a system.

#include <iostream>
#include <mathrix/mathrix.hpp>

using namespace mathrix;

int main() {
    // Quaternion attitude: rotate the body x axis 90 deg about z.
    const Quat q = Quat::fromAxisAngle(Vec3::unitZ(), deg2rad(90.0));
    std::cout << "rotated x axis: " << q.rotate(Vec3::unitX()) << "\n";

    // Monotone interpolation of a thrust curve and its exact integral (total impulse).
    const Interpolator1D thrust({0.0, 0.1, 1.0, 2.0, 2.2}, {0.0, 900.0, 800.0, 600.0, 0.0}, Interp::Pchip,
                                Extrapolation::Zero);
    std::cout << "thrust at 1.5 s: " << thrust(1.5) << " N, impulse " << thrust.integrate(0.0, 2.2) << " Ns\n";

    // Harmonic oscillator with RK4.
    using S = VecN<2>;
    auto f = [](double, const S& y) {
        S d;
        d[0] = y[1];
        d[1] = -y[0];
        return d;
    };
    S y;
    y[0] = 1.0;
    for (int i = 0; i < 100; ++i) y = fixedStep(Integrator::RK4, f, i * 0.01, y, 0.01);
    std::cout << "x(1) = " << y[0] << " (exact " << std::cos(1.0) << ")\n";

    // Linear system and polynomial calculus.
    const Matrix A{{3, 2}, {1, 2}};
    std::cout << "solution of A x = [7, 5]:\n" << A.solve(Matrix::column({7, 5})) << "\n";
    const Polynomial p{1, 0, -3, 1};
    std::cout << "p(x) = " << p << ", p'(x) = " << p.derivative() << ", root in [0,1]: " << p.rootIn(0.0, 1.0) << "\n";
}
