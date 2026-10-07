#pragma once

#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mathrix {

/// Dense, row-major, dynamically sized matrix of doubles.
class Matrix {
public:
    Matrix() = default;
    Matrix(size_t rows, size_t cols, double fill = 0.0) : r_(rows), c_(cols), d_(rows * cols, fill) {}
    Matrix(std::initializer_list<std::initializer_list<double>> rows) {
        r_ = rows.size();
        c_ = r_ ? rows.begin()->size() : 0;
        d_.reserve(r_ * c_);
        for (const auto& row : rows) {
            if (row.size() != c_) throw std::invalid_argument("Matrix: rows must have equal length");
            d_.insert(d_.end(), row.begin(), row.end());
        }
    }

    static Matrix identity(size_t n) {
        Matrix m(n, n);
        for (size_t i = 0; i < n; ++i) m(i, i) = 1.0;
        return m;
    }
    /// Column vector from values.
    static Matrix column(std::initializer_list<double> v) {
        Matrix m(v.size(), 1);
        size_t i = 0;
        for (double x : v) m(i++, 0) = x;
        return m;
    }

    size_t rows() const { return r_; }
    size_t cols() const { return c_; }
    bool isSquare() const { return r_ == c_; }

    double& operator()(size_t r, size_t c) { return d_[r * c_ + c]; }
    double operator()(size_t r, size_t c) const { return d_[r * c_ + c]; }
    double& at(size_t r, size_t c) {
        if (r >= r_ || c >= c_) throw std::out_of_range("Matrix index out of range");
        return (*this)(r, c);
    }
    double at(size_t r, size_t c) const {
        if (r >= r_ || c >= c_) throw std::out_of_range("Matrix index out of range");
        return (*this)(r, c);
    }

    Matrix operator+(const Matrix& o) const { return combine(o, 1.0); }
    Matrix operator-(const Matrix& o) const { return combine(o, -1.0); }
    Matrix operator*(double s) const {
        Matrix m = *this;
        for (double& x : m.d_) x *= s;
        return m;
    }
    friend Matrix operator*(double s, const Matrix& m) { return m * s; }
    Matrix operator*(const Matrix& o) const {
        if (c_ != o.r_) throw std::invalid_argument("Matrix product: dimension mismatch");
        Matrix m(r_, o.c_);
        for (size_t i = 0; i < r_; ++i)
            for (size_t k = 0; k < c_; ++k) {
                const double a = (*this)(i, k);
                if (a == 0.0) continue;
                for (size_t j = 0; j < o.c_; ++j) m(i, j) += a * o(k, j);
            }
        return m;
    }
    bool operator==(const Matrix& o) const { return r_ == o.r_ && c_ == o.c_ && d_ == o.d_; }
    bool operator!=(const Matrix& o) const { return !(*this == o); }

    Matrix transposed() const {
        Matrix m(c_, r_);
        for (size_t i = 0; i < r_; ++i)
            for (size_t j = 0; j < c_; ++j) m(j, i) = (*this)(i, j);
        return m;
    }

    double trace() const {
        requireSquare("trace");
        double t = 0.0;
        for (size_t i = 0; i < r_; ++i) t += (*this)(i, i);
        return t;
    }

    /// Determinant by LU decomposition with partial pivoting.
    double determinant() const {
        requireSquare("determinant");
        Matrix a = *this;
        double det = 1.0;
        for (size_t k = 0; k < r_; ++k) {
            const size_t p = pivot(a, k);
            if (std::abs(a(p, k)) < 1e-300) return 0.0;
            if (p != k) {
                a.swapRows(p, k);
                det = -det;
            }
            det *= a(k, k);
            for (size_t i = k + 1; i < r_; ++i) {
                const double f = a(i, k) / a(k, k);
                for (size_t j = k; j < c_; ++j) a(i, j) -= f * a(k, j);
            }
        }
        return det;
    }

    /// Solves A x = B (B may have several columns) by Gaussian elimination with partial pivoting.
    Matrix solve(const Matrix& b) const {
        requireSquare("solve");
        if (b.r_ != r_) throw std::invalid_argument("Matrix solve: right-hand side has wrong size");
        Matrix a = *this, x = b;
        for (size_t k = 0; k < r_; ++k) {
            const size_t p = pivot(a, k);
            if (std::abs(a(p, k)) < 1e-14 * norm()) throw std::domain_error("Matrix is singular");
            a.swapRows(p, k);
            x.swapRows(p, k);
            for (size_t i = 0; i < r_; ++i) {
                if (i == k) continue;
                const double f = a(i, k) / a(k, k);
                if (f == 0.0) continue;
                for (size_t j = k; j < c_; ++j) a(i, j) -= f * a(k, j);
                for (size_t j = 0; j < x.c_; ++j) x(i, j) -= f * x(k, j);
            }
        }
        for (size_t i = 0; i < r_; ++i)
            for (size_t j = 0; j < x.c_; ++j) x(i, j) /= a(i, i);
        return x;
    }

    /// Inverse (Gauss-Jordan); throws std::domain_error for singular matrices.
    Matrix inverse() const { return solve(identity(r_)); }

    /// Frobenius norm.
    double norm() const {
        double s = 0.0;
        for (double x : d_) s += x * x;
        return std::sqrt(s);
    }

    void swapRows(size_t a, size_t b) {
        if (a == b) return;
        for (size_t j = 0; j < c_; ++j) std::swap((*this)(a, j), (*this)(b, j));
    }

    const std::vector<double>& data() const { return d_; }

private:
    size_t r_ = 0, c_ = 0;
    std::vector<double> d_;

    Matrix combine(const Matrix& o, double s) const {
        if (r_ != o.r_ || c_ != o.c_) throw std::invalid_argument("Matrix sum: dimension mismatch");
        Matrix m = *this;
        for (size_t i = 0; i < d_.size(); ++i) m.d_[i] += s * o.d_[i];
        return m;
    }
    void requireSquare(const char* what) const {
        if (!isSquare()) throw std::invalid_argument(std::string("Matrix ") + what + ": matrix must be square");
    }
    static size_t pivot(const Matrix& a, size_t k) {
        size_t p = k;
        for (size_t i = k + 1; i < a.r_; ++i)
            if (std::abs(a(i, k)) > std::abs(a(p, k))) p = i;
        return p;
    }
};

inline std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (size_t i = 0; i < m.rows(); ++i) {
        os << (i ? " [" : "[[");
        for (size_t j = 0; j < m.cols(); ++j) os << (j ? ", " : "") << m(i, j);
        os << (i + 1 < m.rows() ? "]\n" : "]]");
    }
    return os;
}

}  // namespace mathrix
