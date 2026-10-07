#include <mathrix/mathrix.hpp>

#include "TestFramework.hpp"

using namespace mathrix;

TEST_CASE("quaternion rotation and derivative") {
    const Quat q = Quat::fromAxisAngle(Vec3::unitZ(), kHalfPi);
    const Vec3 v = q.rotate(Vec3::unitX());
    CHECK_NEAR(v.x, 0.0, 1e-12);
    CHECK_NEAR(v.y, 1.0, 1e-12);
    CHECK_NEAR(q.inverseRotate(v).x, 1.0, 1e-12);
    const Quat r = Quat::fromTwoVectors(Vec3::unitX(), Vec3{0.1, 0.2, 0.97}.normalized());
    const Vec3 d = r.rotate(Vec3::unitX());
    CHECK_NEAR(d.dot(Vec3{0.1, 0.2, 0.97}.normalized()), 1.0, 1e-12);
    // Matrix and quaternion rotations agree.
    const Vec3 w{0.3, -0.2, 0.7};
    const Vec3 a = r.rotate(w), b = r.toMatrix() * w;
    CHECK_NEAR((a - b).norm(), 0.0, 1e-12);
    // Integrating qdot for 1 s at 90 deg/s around z gives a 90 deg rotation.
    Quat qi = Quat::identity();
    const Vec3 omega{0.0, 0.0, kHalfPi};
    for (int i = 0; i < 1000; ++i) {
        const double h = 1e-3;
        const Quat k1 = qi.derivative(omega);
        const Quat k2 = (qi + k1 * (h / 2)).derivative(omega);
        const Quat k3 = (qi + k2 * (h / 2)).derivative(omega);
        const Quat k4 = (qi + k3 * h).derivative(omega);
        qi = (qi + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0)).normalized();
    }
    CHECK_NEAR(qi.rotate(Vec3::unitX()).y, 1.0, 1e-9);
}

TEST_CASE("matrix inverse") {
    const Mat3 m{4, 1, 0, 1, 3, 1, 0, 1, 2};
    const Mat3 p = m * m.inverse();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) CHECK_NEAR(p(i, j), i == j ? 1.0 : 0.0, 1e-12);
}

