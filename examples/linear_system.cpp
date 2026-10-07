// Solving a system of n linear equations with n unknowns: A x = b.
//
//    2x +  y -  z =   8
//   -3x -  y + 2z = -11
//   -2x +  y + 2z =  -3        ->  x = 2, y = 3, z = -1
//
// Any size works: build an n x n matrix A and an n x 1 column b, then call A.solve(b).

#include <iostream>
#include <mathrix/mathrix.hpp>

using mathrix::Matrix;

int main() {
    // 1) A small system written out by hand.
    const Matrix A{{2, 1, -1},
                   {-3, -1, 2},
                   {-2, 1, 2}};
    const Matrix b = Matrix::column({8, -11, -3});

    const Matrix x = A.solve(b);
    std::cout << "x =\n" << x << "\n";
    std::cout << "check A*x =\n" << A * x << "\n\n";

    // 2) The same call for n unknowns: here n = 6, A(i,j) = 1 / (i + j + 1) + (i == j ? n : 0).
    const size_t n = 6;
    Matrix An(n, n), bn(n, 1);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) An(i, j) = 1.0 / static_cast<double>(i + j + 1) + (i == j ? double(n) : 0.0);
        bn(i, 0) = static_cast<double>(i + 1);
    }
    const Matrix xn = An.solve(bn);
    std::cout << n << " unknowns, solution =\n" << xn.transposed() << "\n";
    std::cout << "residual |A x - b| = " << (An * xn - bn).norm() << "\n\n";

    // 3) A singular system has no unique solution: solve() throws.
    try {
        Matrix{{1, 2}, {2, 4}}.solve(Matrix::column({3, 6}));
    } catch (const std::exception& e) {
        std::cout << "singular system: " << e.what() << "\n";
    }
}
