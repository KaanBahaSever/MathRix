#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace mathrix {

/// Fixed-size state vector with the arithmetic needed by the ODE integrators.
template <size_t N>
struct VecN {
    std::array<double, N> v{};

    static constexpr size_t size() { return N; }
    double& operator[](size_t i) { return v[i]; }
    double operator[](size_t i) const { return v[i]; }

    VecN operator+(const VecN& o) const { VecN r; for (size_t i = 0; i < N; ++i) r.v[i] = v[i] + o.v[i]; return r; }
    VecN operator-(const VecN& o) const { VecN r; for (size_t i = 0; i < N; ++i) r.v[i] = v[i] - o.v[i]; return r; }
    VecN operator*(double s) const { VecN r; for (size_t i = 0; i < N; ++i) r.v[i] = v[i] * s; return r; }
    VecN& operator+=(const VecN& o) { for (size_t i = 0; i < N; ++i) v[i] += o.v[i]; return *this; }
    friend VecN operator*(double s, const VecN& a) { return a * s; }
};

/// Available integration schemes.
///  - Euler:         1st order explicit (what most first simulators use; needs tiny steps)
///  - Midpoint:      2nd order (RK2 midpoint)
///  - Heun:          2nd order (improved Euler / trapezoidal predictor-corrector)
///  - RK4:           classic 4th order Runge-Kutta (default)
///  - DormandPrince: adaptive 5(4) embedded Runge-Kutta (RK45) with error control
enum class Integrator { Euler, Midpoint, Heun, RK4, DormandPrince };

inline const char* toString(Integrator i) {
    switch (i) {
        case Integrator::Euler: return "euler";
        case Integrator::Midpoint: return "midpoint";
        case Integrator::Heun: return "heun";
        case Integrator::RK4: return "rk4";
        case Integrator::DormandPrince: return "rk45";
    }
    return "rk4";
}

inline Integrator integratorFromString(const std::string& s) {
    if (s == "euler") return Integrator::Euler;
    if (s == "midpoint" || s == "rk2") return Integrator::Midpoint;
    if (s == "heun") return Integrator::Heun;
    if (s == "rk4" || s == "runge-kutta") return Integrator::RK4;
    if (s == "rk45" || s == "dopri" || s == "dormand-prince" || s == "adaptive") return Integrator::DormandPrince;
    throw std::invalid_argument("unknown integrator '" + s + "'");
}

inline bool isAdaptive(Integrator i) { return i == Integrator::DormandPrince; }

/// One fixed step of a single-step scheme. `f(t, y)` returns dy/dt.
template <class State, class F>
State fixedStep(Integrator method, F&& f, double t, const State& y, double h) {
    switch (method) {
        case Integrator::Euler:
            return y + f(t, y) * h;
        case Integrator::Midpoint: {
            const State k1 = f(t, y);
            const State k2 = f(t + 0.5 * h, y + k1 * (0.5 * h));
            return y + k2 * h;
        }
        case Integrator::Heun: {
            const State k1 = f(t, y);
            const State k2 = f(t + h, y + k1 * h);
            return y + (k1 + k2) * (0.5 * h);
        }
        case Integrator::RK4:
        case Integrator::DormandPrince: {
            const State k1 = f(t, y);
            const State k2 = f(t + 0.5 * h, y + k1 * (0.5 * h));
            const State k3 = f(t + 0.5 * h, y + k2 * (0.5 * h));
            const State k4 = f(t + h, y + k3 * h);
            return y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
        }
    }
    return y;
}

/// Result of an embedded Runge-Kutta trial step.
template <class State>
struct EmbeddedStep {
    State y;          ///< 5th order solution
    double error = 0; ///< scaled RMS error (<= 1 means acceptable)
};

/// Dormand-Prince 5(4) trial step with mixed absolute/relative error scaling.
template <class State, class F>
EmbeddedStep<State> dormandPrinceStep(F&& f, double t, const State& y, double h,
                                      double absTol, double relTol) {
    const State k1 = f(t, y);
    const State k2 = f(t + h / 5.0, y + k1 * (h / 5.0));
    const State k3 = f(t + 3.0 * h / 10.0, y + (k1 * (3.0 / 40.0) + k2 * (9.0 / 40.0)) * h);
    const State k4 = f(t + 4.0 * h / 5.0, y + (k1 * (44.0 / 45.0) + k2 * (-56.0 / 15.0) + k3 * (32.0 / 9.0)) * h);
    const State k5 = f(t + 8.0 * h / 9.0,
                       y + (k1 * (19372.0 / 6561.0) + k2 * (-25360.0 / 2187.0) + k3 * (64448.0 / 6561.0) +
                            k4 * (-212.0 / 729.0)) * h);
    const State k6 = f(t + h, y + (k1 * (9017.0 / 3168.0) + k2 * (-355.0 / 33.0) + k3 * (46732.0 / 5247.0) +
                                   k4 * (49.0 / 176.0) + k5 * (-5103.0 / 18656.0)) * h);
    const State y5 = y + (k1 * (35.0 / 384.0) + k3 * (500.0 / 1113.0) + k4 * (125.0 / 192.0) +
                          k5 * (-2187.0 / 6784.0) + k6 * (11.0 / 84.0)) * h;
    const State k7 = f(t + h, y5);
    const State err = (k1 * (35.0 / 384.0 - 5179.0 / 57600.0) + k3 * (500.0 / 1113.0 - 7571.0 / 16695.0) +
                       k4 * (125.0 / 192.0 - 393.0 / 640.0) + k5 * (-2187.0 / 6784.0 + 92097.0 / 339200.0) +
                       k6 * (11.0 / 84.0 - 187.0 / 2100.0) + k7 * (-1.0 / 40.0)) * h;
    double sum = 0.0;
    for (size_t i = 0; i < State::size(); ++i) {
        const double sc = absTol + relTol * std::max(std::abs(y[i]), std::abs(y5[i]));
        const double e = err[i] / sc;
        sum += e * e;
    }
    return {y5, std::sqrt(sum / static_cast<double>(State::size()))};
}

/// Standard step-size controller for a 5th order method.
inline double nextStepSize(double h, double error, double minFactor = 0.2, double maxFactor = 5.0) {
    if (error <= 0.0) return h * maxFactor;
    const double factor = 0.9 * std::pow(error, -0.2);
    return h * std::clamp(factor, minFactor, maxFactor);
}

}  // namespace mathrix
