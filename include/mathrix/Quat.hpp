#pragma once

#include <cmath>

#include "Mat3.hpp"
#include "Vec3.hpp"

namespace mathrix {

/// Unit quaternion (Hamilton convention, w + xi + yj + zk) representing a rotation.
/// `q.rotate(v)` maps a vector expressed in the body frame into the reference frame.
struct Quat {
    double w = 1.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Quat() = default;
    constexpr Quat(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}

    static constexpr Quat identity() { return {1, 0, 0, 0}; }

    static Quat fromAxisAngle(const Vec3& axis, double angle) {
        const Vec3 a = axis.normalized();
        const double s = std::sin(0.5 * angle);
        return {std::cos(0.5 * angle), a.x * s, a.y * s, a.z * s};
    }

    /// Shortest rotation taking unit vector `from` onto unit vector `to`.
    static Quat fromTwoVectors(const Vec3& from, const Vec3& to) {
        const Vec3 f = from.normalized();
        const Vec3 t = to.normalized();
        const double d = f.dot(t);
        if (d < -0.999999999) {
            // Opposite vectors: rotate 180 deg around any orthogonal axis.
            Vec3 axis = Vec3::unitX().cross(f);
            if (axis.norm2() < 1e-12) axis = Vec3::unitY().cross(f);
            return fromAxisAngle(axis, 3.14159265358979323846);
        }
        const Vec3 c = f.cross(t);
        Quat q{1.0 + d, c.x, c.y, c.z};
        return q.normalized();
    }

    constexpr Quat operator*(const Quat& o) const {
        return {w * o.w - x * o.x - y * o.y - z * o.z,
                w * o.x + x * o.w + y * o.z - z * o.y,
                w * o.y - x * o.z + y * o.w + z * o.x,
                w * o.z + x * o.y - y * o.x + z * o.w};
    }
    constexpr Quat operator+(const Quat& o) const { return {w + o.w, x + o.x, y + o.y, z + o.z}; }
    constexpr Quat operator*(double s) const { return {w * s, x * s, y * s, z * s}; }

    constexpr Quat conjugate() const { return {w, -x, -y, -z}; }
    double norm() const { return std::sqrt(w * w + x * x + y * y + z * z); }
    Quat normalized() const {
        const double n = norm();
        return n > 1e-300 ? Quat{w / n, x / n, y / n, z / n} : identity();
    }

    /// Rotate body-frame vector into reference frame: v' = q v q*.
    constexpr Vec3 rotate(const Vec3& v) const {
        // Optimised form of q * (0,v) * q^-1 for unit quaternions.
        const Vec3 u{x, y, z};
        const Vec3 t = u.cross(v) * 2.0;
        return v + t * w + u.cross(t);
    }
    /// Rotate reference-frame vector into body frame.
    constexpr Vec3 inverseRotate(const Vec3& v) const { return conjugate().rotate(v); }

    /// Time derivative for body-frame angular velocity `omegaBody`: qdot = 0.5 q (0, w).
    constexpr Quat derivative(const Vec3& omegaBody) const {
        return Quat{w, x, y, z} * Quat{0.0, omegaBody.x, omegaBody.y, omegaBody.z} * 0.5;
    }

    /// Rotation matrix (body -> reference).
    constexpr Mat3 toMatrix() const {
        return {1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
                2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
                2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)};
    }
};

}  // namespace mathrix
