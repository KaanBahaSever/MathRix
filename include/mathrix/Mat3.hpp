#pragma once

#include <array>
#include <cmath>

#include "Vec3.hpp"

namespace mathrix {

/// Row-major 3x3 matrix.
struct Mat3 {
    std::array<double, 9> m{0, 0, 0, 0, 0, 0, 0, 0, 0};

    constexpr Mat3() = default;
    constexpr Mat3(double a00, double a01, double a02,
                   double a10, double a11, double a12,
                   double a20, double a21, double a22)
        : m{a00, a01, a02, a10, a11, a12, a20, a21, a22} {}

    static constexpr Mat3 identity() { return {1, 0, 0, 0, 1, 0, 0, 0, 1}; }
    static constexpr Mat3 diag(double a, double b, double c) { return {a, 0, 0, 0, b, 0, 0, 0, c}; }
    static constexpr Mat3 fromRows(const Vec3& r0, const Vec3& r1, const Vec3& r2) {
        return {r0.x, r0.y, r0.z, r1.x, r1.y, r1.z, r2.x, r2.y, r2.z};
    }
    static constexpr Mat3 fromColumns(const Vec3& c0, const Vec3& c1, const Vec3& c2) {
        return {c0.x, c1.x, c2.x, c0.y, c1.y, c2.y, c0.z, c1.z, c2.z};
    }

    constexpr double operator()(int r, int c) const { return m[static_cast<size_t>(r * 3 + c)]; }
    double& operator()(int r, int c) { return m[static_cast<size_t>(r * 3 + c)]; }

    constexpr Vec3 row(int r) const { return {(*this)(r, 0), (*this)(r, 1), (*this)(r, 2)}; }
    constexpr Vec3 col(int c) const { return {(*this)(0, c), (*this)(1, c), (*this)(2, c)}; }

    constexpr Vec3 operator*(const Vec3& v) const {
        return {m[0] * v.x + m[1] * v.y + m[2] * v.z,
                m[3] * v.x + m[4] * v.y + m[5] * v.z,
                m[6] * v.x + m[7] * v.y + m[8] * v.z};
    }

    constexpr Mat3 operator*(const Mat3& o) const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0;
                for (int k = 0; k < 3; ++k) s += (*this)(i, k) * o(k, j);
                r.m[static_cast<size_t>(i * 3 + j)] = s;
            }
        return r;
    }

    constexpr Mat3 operator*(double s) const {
        Mat3 r;
        for (size_t i = 0; i < 9; ++i) r.m[i] = m[i] * s;
        return r;
    }
    constexpr Mat3 operator+(const Mat3& o) const {
        Mat3 r;
        for (size_t i = 0; i < 9; ++i) r.m[i] = m[i] + o.m[i];
        return r;
    }
    constexpr Mat3 operator-(const Mat3& o) const {
        Mat3 r;
        for (size_t i = 0; i < 9; ++i) r.m[i] = m[i] - o.m[i];
        return r;
    }

    constexpr Mat3 transposed() const {
        return {m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]};
    }

    constexpr double determinant() const {
        return m[0] * (m[4] * m[8] - m[5] * m[7]) -
               m[1] * (m[3] * m[8] - m[5] * m[6]) +
               m[2] * (m[3] * m[7] - m[4] * m[6]);
    }

    /// Inverse; returns zero matrix when singular.
    Mat3 inverse() const {
        const double det = determinant();
        if (std::abs(det) < 1e-300) return Mat3{};
        const double id = 1.0 / det;
        return Mat3{(m[4] * m[8] - m[5] * m[7]) * id, (m[2] * m[7] - m[1] * m[8]) * id, (m[1] * m[5] - m[2] * m[4]) * id,
                    (m[5] * m[6] - m[3] * m[8]) * id, (m[0] * m[8] - m[2] * m[6]) * id, (m[2] * m[3] - m[0] * m[5]) * id,
                    (m[3] * m[7] - m[4] * m[6]) * id, (m[1] * m[6] - m[0] * m[7]) * id, (m[0] * m[4] - m[1] * m[3]) * id};
    }
};

}  // namespace mathrix
