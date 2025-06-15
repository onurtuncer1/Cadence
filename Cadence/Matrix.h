//============================================================================//
// Cadence - Constexpr Math for CAD Kernels
// SPDX-License-Identifier: MIT
// Author: Onur Tuncer, PhD <onur.tuncer@itu.edu.tr>
// File: [filename]
// Desc: [1-2 line description of purpose]
//============================================================================//


#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <cmath>
#include "Vector.h"

namespace Cadence {

template<std::floating_point T, std::size_t Rows, std::size_t Cols>
class Matrix {
    std::array<T, Rows * Cols> data_{};

public:
    static constexpr std::size_t row_count = Rows;
    static constexpr std::size_t col_count = Cols;

    constexpr Matrix() = default;

    constexpr Matrix(const std::array<T, Rows * Cols>& values)
        : data_{values} {}

    constexpr T& operator()(std::size_t row, std::size_t col) {
        return data_[row * Cols + col];
    }

    constexpr const T& operator()(std::size_t row, std::size_t col) const {
        return data_[row * Cols + col];
    }

    constexpr std::size_t rows() const { return Rows; }
    constexpr std::size_t cols() const { return Cols; }

    constexpr const std::array<T, Rows * Cols>& raw() const { return data_; }

    // Identity matrix (only for square matrices)
    static constexpr Matrix identity() requires (Rows == Cols) {
        Matrix m;
        for (std::size_t i = 0; i < Rows; ++i)
            m(i, i) = T{1};
        return m;
    }

    // Matrix-matrix multiplication
    template<std::size_t OtherCols>
    constexpr Matrix<T, Rows, OtherCols> operator*(const Matrix<T, Cols, OtherCols>& rhs) const {
        Matrix<T, Rows, OtherCols> result;
        for (std::size_t i = 0; i < Rows; ++i) {
            for (std::size_t j = 0; j < OtherCols; ++j) {
                T sum = T{};
                for (std::size_t k = 0; k < Cols; ++k) {
                    sum += (*this)(i, k) * rhs(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    // Matrix-vector multiplication
    constexpr Vector<T, Rows> operator*(const Vector<T, Cols>& vec) const {
        Vector<T, Rows> result;
        for (std::size_t i = 0; i < Rows; ++i) {
            T sum = T{};
            for (std::size_t j = 0; j < Cols; ++j) {
                sum += (*this)(i, j) * vec[j];
            }
            result[i] = sum;
        }
        return result;
    }

    // Scalar ops
    constexpr Matrix operator*(T scalar) const {
        Matrix result;
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            result.data_[i] = data_[i] * scalar;
        return result;
    }

    constexpr Matrix operator/(T scalar) const {
        Matrix result;
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            result.data_[i] = data_[i] / scalar;
        return result;
    }

    constexpr Matrix operator+(const Matrix& other) const {
        Matrix result;
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            result.data_[i] = data_[i] + other.data_[i];
        return result;
    }

    constexpr Matrix operator-(const Matrix& other) const {
        Matrix result;
        for (std::size_t i = 0; i < Rows * Cols; ++i)
            result.data_[i] = data_[i] - other.data_[i];
        return result;
    }

    constexpr T determinant() const requires (Rows == Cols) {
    return Cadence::determinant(*this);
}

constexpr bool has_nan() const {
    for (std::size_t i = 0; i < Rows * Cols; ++i) {
        if (std::isnan(data_[i])) return true;
    }
    return false;
}


constexpr T max_element() const {
    T max = data_[0];
    for (std::size_t i = 1; i < Rows * Cols; ++i) {
        if (data_[i] > max)
            max = data_[i];
    }
    return max;
}



};

// ==========================================
// Eigen-style compile-time insertion support
// ==========================================

template<std::floating_point T, std::size_t Rows, std::size_t Cols>
struct MatrixInserter {
    std::array<T, Rows * Cols> buffer{};
    std::size_t index = 0;

    constexpr MatrixInserter& operator,(T value) {
        buffer[index++] = value;
        return *this;
    }

    constexpr operator Matrix<T, Rows, Cols>() const {
        Matrix<T, Rows, Cols> mat;
        for (std::size_t i = 0; i < Rows; ++i)
            for (std::size_t j = 0; j < Cols; ++j)
                mat(i, j) = buffer[i * Cols + j];
        return mat;
    }
};

template<std::floating_point T, std::size_t Rows, std::size_t Cols>
constexpr MatrixInserter<T, Rows, Cols> operator<<(const Matrix<T, Rows, Cols>&, T value) {
    MatrixInserter<T, Rows, Cols> inserter{};
    inserter.buffer[0] = value;
    inserter.index = 1;
    return inserter;
}



template<std::floating_point T, std::size_t N>
constexpr T determinant(const Matrix<T, N, N>& m);

// Base case: 1x1
template<std::floating_point T>
constexpr T determinant(const Matrix<T, 1, 1>& m) {
    return m(0, 0);
}

// Base case: 2x2
template<std::floating_point T>
constexpr T determinant(const Matrix<T, 2, 2>& m) {
    return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
}

// Base case: 3x3 (hardcoded for performance)
template<std::floating_point T>
constexpr T determinant(const Matrix<T, 3, 3>& m) {
    return
        m(0,0) * (m(1,1) * m(2,2) - m(1,2) * m(2,1)) -
        m(0,1) * (m(1,0) * m(2,2) - m(1,2) * m(2,0)) +
        m(0,2) * (m(1,0) * m(2,1) - m(1,1) * m(2,0));
}

// Helper to compute the minor matrix (submatrix after removing row r and col c)
template<std::floating_point T, std::size_t N>
constexpr Matrix<T, N - 1, N - 1> minor_matrix(const Matrix<T, N, N>& m, std::size_t row, std::size_t col) {
    Matrix<T, N - 1, N - 1> result;
    for (std::size_t i = 0, ri = 0; i < N; ++i) {
        if (i == row) continue;
        for (std::size_t j = 0, rj = 0; j < N; ++j) {
            if (j == col) continue;
            result(ri, rj++) = m(i, j);
        }
        ++ri;
    }
    return result;
}

// Recursive case for N > 3
template<std::floating_point T, std::size_t N>
constexpr T determinant(const Matrix<T, N, N>& m) {
    T det = T{0};
    for (std::size_t col = 0; col < N; ++col) {
        const T sign = (col % 2 == 0) ? T{1} : T{-1};
        det += sign * m(0, col) * determinant(minor_matrix(m, 0, col));
    }
    return det;
}

template<typename T, std::size_t N>
constexpr bool is_symmetric(const Matrix<T, N, N>& A, T tol = 1e-9) {
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j)
            if (std::abs(A(i, j) - A(j, i)) > tol)
                return false;
    return true;
}

template<typename T, std::size_t N>
constexpr bool is_positive_definite(const Matrix<T, N, N>& A, T tol = 1e-9) {
    for (std::size_t i = 0; i < N; ++i) {
        T sum_diag = 0;
        for (std::size_t k = 0; k < i; ++k)
            sum_diag += A(i, k) * A(i, k);
        T val = A(i, i) - sum_diag;
        if (val <= tol)
            return false;
    }
    return true;
}

// Concept with matrix value passed in
template<typename T, std::size_t N>
concept SymmetricPositiveDefinite = requires(const Matrix<T, N, N>& A) {
    { is_symmetric(A) } -> std::same_as<bool>;
    { is_positive_definite(A) } -> std::same_as<bool>;
};

template<std::floating_point T, std::size_t N>
requires SymmetricPositiveDefinite<T, N>
constexpr Matrix<T, N, N> cholesky(const Matrix<T, N, N>& A) {
    Matrix<T, N, N> L;

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            T sum = 0;
            for (std::size_t k = 0; k < j; ++k)
                sum += L(i, k) * L(j, k);

            if (i == j) {
                T val = A(i, i) - sum;
                L(i, j) = std::sqrt(val);
            } else {
                L(i, j) = (A(i, j) - sum) / L(j, j);
            }
        }
    }

    return L;
}

// Construct a diagonal matrix from a Vector<T, N>
static constexpr Matrix from_diagonal(const Vector<T, (Rows < Cols ? Rows : Cols)>& diag) {
    Matrix result{};
    for (std::size_t i = 0; i < diag.size(); ++i)
        result(i, i) = diag[i];
    return result;
}

// Transpose of a matrix
constexpr Matrix<T, Cols, Rows> transpose() const {
    Matrix<T, Cols, Rows> result;
    for (std::size_t i = 0; i < Rows; ++i)
        for (std::size_t j = 0; j < Cols; ++j)
            result(j, i) = (*this)(i, j);
    return result;
}

// Trace (sum of diagonal elements)
constexpr T trace() const requires (Rows == Cols) {
    T sum = T{0};
    for (std::size_t i = 0; i < Rows; ++i)
        sum += (*this)(i, i);
    return sum;
}

// All-zero matrix
static constexpr Matrix zero() {
    return Matrix{}; // default-initialized to 0
}

// Fill with a constant value
static constexpr Matrix filled(T value) {
    Matrix result;
    for (std::size_t i = 0; i < Rows * Cols; ++i)
        result.data_[i] = value;
    return result;
}

// Trait to detect a matrix expression (future-proof)
template<typename T>
concept MatrixExpression = requires(T expr, std::size_t i, std::size_t j) {
    { expr(i, j) } -> std::convertible_to<double>;
    { expr.rows() } -> std::convertible_to<std::size_t>;
    { expr.cols() } -> std::convertible_to<std::size_t>;
};


// ==========================================
// Common type aliases
// ==========================================

template<typename T> using Matrix2x2 = Matrix<T, 2, 2>;
template<typename T> using Matrix3x3 = Matrix<T, 3, 3>;
template<typename T> using Matrix4x4 = Matrix<T, 4, 4>;

} // namespace Cadence
