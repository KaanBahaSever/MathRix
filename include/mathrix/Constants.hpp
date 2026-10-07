#pragma once

namespace mathrix {

inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kTwoPi = 2.0 * kPi;
inline constexpr double kHalfPi = 0.5 * kPi;
inline constexpr double kDegToRad = kPi / 180.0;
inline constexpr double kRadToDeg = 180.0 / kPi;

inline constexpr double deg2rad(double deg) { return deg * kDegToRad; }
inline constexpr double rad2deg(double rad) { return rad * kRadToDeg; }

template <class T>
constexpr T clamp(T v, T lo, T hi) { return v < lo ? lo : (v > hi ? hi : v); }

template <class T>
constexpr T sq(T v) { return v * v; }

template <class T>
constexpr int sign(T v) { return (T(0) < v) - (v < T(0)); }

/// Smooth Hermite step: 0 for x <= e0, 1 for x >= e1, C1-continuous in between.
inline double smoothstep(double e0, double e1, double x) {
    if (x <= e0) return 0.0;
    if (x >= e1) return 1.0;
    const double t = (x - e0) / (e1 - e0);
    return t * t * (3.0 - 2.0 * t);
}

inline double lerp(double a, double b, double t) { return a + (b - a) * t; }

}  // namespace mathrix
