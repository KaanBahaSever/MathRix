#pragma once

#include <cmath>
#include <ostream>

namespace mathrix {

/// Plain 3-component double vector. Trivially copyable, no hidden allocations.
struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    static constexpr Vec3 zero() { return {0.0, 0.0, 0.0}; }
    static constexpr Vec3 unitX() { return {1.0, 0.0, 0.0}; }
    static constexpr Vec3 unitY() { return {0.0, 1.0, 0.0}; }
    static constexpr Vec3 unitZ() { return {0.0, 0.0, 1.0}; }

    constexpr double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
    double& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }

    constexpr Vec3 operator-() const { return {-x, -y, -z}; }
    constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
    Vec3& operator/=(double s) { x /= s; y /= s; z /= s; return *this; }

    constexpr bool operator==(const Vec3& o) const { return x == o.x && y == o.y && z == o.z; }
    constexpr bool operator!=(const Vec3& o) const { return !(*this == o); }

    constexpr double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    constexpr Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    constexpr double norm2() const { return x * x + y * y + z * z; }
    double norm() const { return std::sqrt(norm2()); }

    /// Unit vector in the same direction; returns zero for (near) zero vectors.
    Vec3 normalized() const {
        const double n = norm();
        return n > 1e-300 ? (*this) / n : Vec3{};
    }

    bool isFinite() const { return std::isfinite(x) && std::isfinite(y) && std::isfinite(z); }
};

constexpr Vec3 operator*(double s, const Vec3& v) { return v * s; }
constexpr double dot(const Vec3& a, const Vec3& b) { return a.dot(b); }
constexpr Vec3 cross(const Vec3& a, const Vec3& b) { return a.cross(b); }
inline double norm(const Vec3& v) { return v.norm(); }

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
}

}  // namespace mathrix
