#include <mathrix/mathrix.hpp>

#include "TestFramework.hpp"

using namespace mathrix;

TEST_CASE("matrix: arithmetic, transpose and product") {
    const Matrix a{{1, 2, 3}, {4, 5, 6}};
    const Matrix b{{7, 8}, {9, 10}, {11, 12}};
    const Matrix c = a * b;
    CHECK(c.rows() == 2 && c.cols() == 2);
    CHECK(c == (Matrix{{58, 64}, {139, 154}}));
    CHECK(a.transposed() == (Matrix{{1, 4}, {2, 5}, {3, 6}}));
    CHECK((a + a) == a * 2.0);
    CHECK((a - a).norm() == 0.0);
    CHECK_THROWS(a * a);
    CHECK_THROWS(a + b);
    CHECK_THROWS(a.at(2, 0));
}

TEST_CASE("matrix: determinant, inverse and solve") {
    const Matrix m{{4, -2, 1}, {-2, 4, -2}, {1, -2, 4}};
    CHECK_NEAR(m.determinant(), 36.0, 1e-12);
    CHECK_NEAR(m.trace(), 12.0, 1e-12);
    const Matrix p = m * m.inverse();
    CHECK_NEAR((p - Matrix::identity(3)).norm(), 0.0, 1e-12);
    const Matrix x = m.solve(Matrix::column({11, -16, 17}));
    CHECK_NEAR(x(0, 0), 1.0, 1e-12);
    CHECK_NEAR(x(1, 0), -2.0, 1e-12);
    CHECK_NEAR(x(2, 0), 3.0, 1e-12);
    // Pivoting is required here (zero on the diagonal).
    const Matrix q{{0, 1}, {1, 0}};
    CHECK_NEAR(q.determinant(), -1.0, 1e-12);
    CHECK(q.inverse() == q);
    const Matrix singular{{1, 2}, {2, 4}};
    CHECK_NEAR(singular.determinant(), 0.0, 1e-12);
    CHECK_THROWS(singular.inverse());
}

TEST_CASE("polynomial: evaluation, calculus and algebra") {
    const Polynomial p{1, -3, 0, 2};  // 2x^3 - 3x + 1
    CHECK(p.degree() == 3);
    CHECK_NEAR(p(2.0), 11.0, 1e-12);
    CHECK(p.derivative() == (Polynomial{-3, 0, 6}));
    CHECK_NEAR(p.integrate(0.0, 2.0), 2.0 * 16.0 / 4.0 - 3.0 * 4.0 / 2.0 + 2.0, 1e-12);
    CHECK(p.integral(5.0).derivative() == p);
    const Polynomial q{-1, 1};  // x - 1
    CHECK((p * q)(3.0) == p(3.0) * q(3.0));
    CHECK((p + q) == (Polynomial{0, -2, 0, 2}));
    CHECK((p - p).degree() == -1);
    // x^2 - 2 has a root at sqrt(2).
    const Polynomial r{-2, 0, 1};
    CHECK_NEAR(r.rootIn(0.0, 2.0), std::sqrt(2.0), 1e-12);
}
