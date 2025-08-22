//------------------------------------------------------------------------------
// Project: Cadence Math Library
// Copyright (c) 2025, Onur Tuncer, PhD, Istanbul Technical University
// SPDX-License-Identifier: MIT
// License-Filename: LICENSE
//------------------------------------------------------------------------------

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <cmath>
#include "Vector.h"

/**
 * @file Matrix.h
 * @brief Provides a template class for matrix operations.
 *
 * This file defines a Matrix class that supports various operations such as
 * addition, subtraction, multiplication, and determinant calculation. It also
 * includes functionality for creating identity matrices, checking for NaN values,
 * and finding the maximum element in the matrix.
 */

namespace Cadence
{

// Forward-declare Matrix so we can name it in the function prototype
template<std::floating_point T, std::size_t Rows, std::size_t Cols>
class Matrix;

// Forward-declare the free-function determinant used by Matrix::determinant()
template<std::floating_point T, std::size_t N>
constexpr T determinant(const Matrix<T, N, N>& m);

struct NoInit{};

/**
 * @class Matrix
 * @brief A class representing a mathematical matrix.
 *
 * @tparam T The type of the matrix elements, which must be a floating-point type.
 * @tparam Rows The number of rows in the matrix.
 * @tparam Cols The number of columns in the matrix.
 *
 * The Matrix class provides methods for matrix arithmetic, including:
 * - Element access via operator()
 * - Matrix-matrix multiplication
 * - Matrix-vector multiplication
 * - Scalar multiplication and division
 * - Addition and subtraction of matrices
 * - Determinant calculation for square matrices
 * - Identity matrix generation
 * - Checking for NaN values
 * - Finding the maximum element
 */
template<std::floating_point T, std::size_t Rows, std::size_t Cols>
class Matrix
{
	std::array<T, Rows *Cols> data_{};

public:
	static constexpr std::size_t row_count = Rows;
	static constexpr std::size_t col_count = Cols;

	/** @brief Default constructor */
	constexpr Matrix() = default;

    /** @brief Constructor without initialization */
	constexpr Matrix(NoInit) noexcept {}

	/** @brief Construct from array of values */
	constexpr Matrix(const std::array<T, Rows *Cols> &values)
		: data_{values} {}

	/** @brief Element access */
	constexpr T &operator()(std::size_t row, std::size_t col)
	{
		return data_[row * Cols + col];
	}

	/** @brief Const element access */
	constexpr const T &operator()(std::size_t row, std::size_t col) const
	{
		return data_[row * Cols + col];
	}

	/** @brief Get number of rows */
	constexpr std::size_t rows() const { return Rows; }

	/** @brief Get number of columns */
	constexpr std::size_t cols() const { return Cols; }

	/** @brief Access raw data */
	constexpr const std::array<T, Rows *Cols> &raw() const { return data_; }

	/** @brief Create identity matrix (square only) */
	static constexpr Matrix identity() requires(Rows == Cols)
	{
		Matrix m;

		for (std::size_t i = 0; i < Rows; ++i)
			m(i, i) = T{1};

		return m;
	}

	/** @brief Matrix multiplication */
	template<std::size_t OtherCols>
	constexpr Matrix<T, Rows, OtherCols> operator*(const Matrix<T, Cols, OtherCols> &rhs) const
	{
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

	/** @brief Matrix-vector multiplication */
	constexpr Vector<T, Rows> operator*(const Vector<T, Cols> &vec) const
	{
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

	/** @brief Scalar multiplication */
	constexpr Matrix operator*(T scalar) const
	{
		Matrix result;

		for (std::size_t i = 0; i < Rows * Cols; ++i) {
			result.data_[i] = data_[i] * scalar;
		}

		return result;
	}

	/** @brief Scalar division */
	constexpr Matrix operator/(T scalar) const
	{
		Matrix result;

		for (std::size_t i = 0; i < Rows * Cols; ++i) {
			result.data_[i] = data_[i] / scalar;
		}

		return result;
	}

	/** @brief Matrix addition */
	constexpr Matrix operator+(const Matrix &other) const
	{
		Matrix result;

		for (std::size_t i = 0; i < Rows * Cols; ++i) {
			result.data_[i] = data_[i] + other.data_[i];
		}

		return result;
	}

	/** @brief Matrix subtraction */
	constexpr Matrix operator-(const Matrix &other) const
	{
		Matrix result;

		for (std::size_t i = 0; i < Rows * Cols; ++i) {
			result.data_[i] = data_[i] - other.data_[i];
		}

		return result;
	}

	/** @brief Determinant (square only) */     //TODO [Onur] looks like we need a forward decleration of determinant 
	constexpr T determinant() const requires(Rows == Cols)
	{
		return Cadence::determinant(*this);
	}

	/** @brief Check for NaNs */
	constexpr bool has_nan() const
	{
		for (std::size_t i = 0; i < Rows * Cols; ++i) {
			if (std::isnan(data_[i])) { return true; }
		}

		return false;
	}

	/** @brief Maximum element */
	constexpr T max_element() const
	{
		T max = data_[0];

		for (std::size_t i = 1; i < Rows * Cols; ++i) {
			if (data_[i] > max) {
				max = data_[i];
			}
		}

		return max;
	}

	/** @brief Transpose matrix */
	constexpr Matrix<T, Cols, Rows> transpose() const
	{
		Matrix<T, Cols, Rows> result;

		for (std::size_t i = 0; i < Rows; ++i)
			for (std::size_t j = 0; j < Cols; ++j) {
				result(j, i) = (*this)(i, j);
			}

		return result;
	}

	/** @brief Trace of matrix (sum of diagonal elements) */
	constexpr T trace() const requires(Rows == Cols)
	{
		T sum = T{0};

		for (std::size_t i = 0; i < Rows; ++i) {
			sum += (*this)(i, i);
		}

		return sum;
	}

	/** @brief Zero matrix */
	static constexpr Matrix zero()
	{
		return Matrix{};
	}

	/** @brief Matrix filled with a constant */
	static constexpr Matrix filled(T value)
	{
		Matrix result;

		for (std::size_t i = 0; i < Rows * Cols; ++i) {
			result.data_[i] = value;
		}

		return result;
	}

	/** @brief Diagonal matrix from vector */
	static constexpr Matrix from_diagonal(const Vector < T, (Rows < Cols ? Rows : Cols) > & diag)
	{
		Matrix result{};

		for (std::size_t i = 0; i < diag.size(); ++i) {
			result(i, i) = diag[i];
		}

		return result;
	}
};

/** @brief Eigen-style compile-time insertion */
template<std::floating_point T, std::size_t Rows, std::size_t Cols>
struct MatrixInserter {
	std::array<T, Rows *Cols> buffer{};
	std::size_t index = 0;

	constexpr MatrixInserter &operator, (T value)
	{
		buffer[index++] = value;
		return *this;
	}

	constexpr operator Matrix<T, Rows, Cols>() const
	{
		Matrix<T, Rows, Cols> mat;

		for (std::size_t i = 0; i < Rows; ++i)
			for (std::size_t j = 0; j < Cols; ++j) {
				mat(i, j) = buffer[i * Cols + j];
			}

		return mat;
	}
};

/** @brief Start matrix insertion */
template<std::floating_point T, std::size_t Rows, std::size_t Cols>
constexpr MatrixInserter<T, Rows, Cols> operator<<(const Matrix<T, Rows, Cols> &, T value)
{
	MatrixInserter<T, Rows, Cols> inserter{};
	inserter.buffer[0] = value;
	inserter.index = 1;
	return inserter;
}

/** @brief Determinant helpers and implementation */
template<std::floating_point T, std::size_t N>
constexpr T determinant(const Matrix<T, N, N> &m);

template<std::floating_point T>
constexpr T determinant(const Matrix<T, 1, 1> &m) { return m(0, 0); }

template<std::floating_point T>
constexpr T determinant(const Matrix<T, 2, 2> &m)
{
	return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
}

template<std::floating_point T>
constexpr T determinant(const Matrix<T, 3, 3> &m)
{
	return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) -
		   m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
		   m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
}

/** @brief Minor matrix helper */
template<std::floating_point T, std::size_t N>
constexpr Matrix < T, N - 1, N - 1 > minor_matrix(const Matrix<T, N, N> &m, std::size_t row, std::size_t col)
{
	Matrix < T, N - 1, N - 1 > result;

	for (std::size_t i = 0, ri = 0; i < N; ++i) {
		if (i == row) { continue; }

		for (std::size_t j = 0, rj = 0; j < N; ++j) {
			if (j == col) { continue; }

			result(ri, rj++) = m(i, j);
		}

		++ri;
	}

	return result;
}

/** @brief General determinant for N > 3 */
template<std::floating_point T, std::size_t N>
constexpr T determinant(const Matrix<T, N, N> &m)
{
	T det = T{0};

	for (std::size_t col = 0; col < N; ++col) {
		const T sign = (col % 2 == 0) ? T{1} : T{-1};
		det += sign * m(0, col) * determinant(minor_matrix(m, 0, col));
	}

	return det;
}

/** @brief Symmetry check */
template<typename T, std::size_t N>
constexpr bool is_symmetric(const Matrix<T, N, N> &A, T tol = 1e-9)
{
	for (std::size_t i = 0; i < N; ++i)
		for (std::size_t j = i + 1; j < N; ++j)
			if (std::abs(A(i, j) - A(j, i)) > tol) {
				return false;
			}

	return true;
}

/** @brief Positive definiteness check */
template<typename T, std::size_t N>
constexpr bool is_positive_definite(const Matrix<T, N, N> &A, T tol = 1e-9)
{
	for (std::size_t i = 0; i < N; ++i) {
		T sum_diag = 0;

		for (std::size_t k = 0; k < i; ++k) {
			sum_diag += A(i, k) * A(i, k);
		}

		T val = A(i, i) - sum_diag;

		if (val <= tol) {
			return false;
		}
	}

	return true;
}

/** @brief Concept for SPD matrix */
template<typename T, std::size_t N>
concept SymmetricPositiveDefinite = requires(const Matrix<T, N, N> &A)
{
	{ is_symmetric(A) } -> std::same_as<bool>;
	{ is_positive_definite(A) } -> std::same_as<bool>;
};

/** @brief Cholesky decomposition */
template<std::floating_point T, std::size_t N>
requires SymmetricPositiveDefinite<T, N>
constexpr Matrix<T, N, N> cholesky(const Matrix<T, N, N> &A)
{
	Matrix<T, N, N> L;

	for (std::size_t i = 0; i < N; ++i) {
		for (std::size_t j = 0; j <= i; ++j) {
			T sum = 0;

			for (std::size_t k = 0; k < j; ++k) {
				sum += L(i, k) * L(j, k);
			}

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

/** @brief Concept to detect a matrix expression */
template<typename T>
concept MatrixExpression = requires(T expr, std::size_t i, std::size_t j)
{
	{ expr(i, j) } -> std::convertible_to<double>;
	{ expr.rows() } -> std::convertible_to<std::size_t>;
	{ expr.cols() } -> std::convertible_to<std::size_t>;
};

/** @brief Type alias for a 2x2 matrix */
 template<typename T> using Matrix2x2 = Matrix<T, 2, 2>;
/** @brief Type alias for a 3x3 matrix */
template<typename T> using Matrix3x3 = Matrix<T, 3, 3>;
/** @brief Type alias for a 4x4 matrix */
template<typename T> using Matrix4x4 = Matrix<T, 4, 4>;

} // namespace Cadence
