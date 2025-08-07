//------------------------------------------------------------------------------
// Project: Cadence Math Library
// Copyright (c) 2025, Onur Tuncer, PhD, Istanbul Technical University
// SPDX-License-Identifier: MIT
// License-Filename: LICENSE
//------------------------------------------------------------------------------

#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>

/**
 * @file Vector.h
 * @brief Provides a template class for fixed-size vector operations.
 */

namespace Cadence
{

/**
 * @class Vector
 * @brief A mathematical vector of fixed size.
 *
 * @tparam T Element type, must satisfy std::floating_point.
 * @tparam N Dimension of the vector.
 */
template<std::floating_point T, std::size_t N>
class Vector
{
	std::array<T, N> data_{};

public:
	/// Static size of the vector
	static constexpr std::size_t static_size = N;

	/** @brief Default constructor */
	constexpr Vector() = default;

	/**
	 * @brief Variadic constructor
	 * @tparam Args Parameter pack for element initialization
	 */
	template<typename... Args>
	requires(sizeof...(Args) == N) &&(std::conjunction_v<std::is_convertible<Args, T>...>)
	constexpr Vector(Args... args) : data_{static_cast<T>(args)...} {}

	/**
	 * @brief Construct from std::array
	 * @param init The array to initialize from
	 */
	constexpr Vector(std::array<T, N> init) : data_{init} {}

	/** @brief Element access (mutable) */
	constexpr T &operator[](std::size_t i) { return data_[i]; }

	/** @brief Element access (const) */
	constexpr const T &operator[](std::size_t i) const { return data_[i]; }

	/** @brief Returns the number of elements */
	constexpr std::size_t size() const { return N; }

	// Named accessors (like Eigen)
	/// @brief Access x component (if N > 0)
	constexpr T &x() requires(N > 0) { return data_[0]; }
	/// @brief Access y component (if N > 1)
	constexpr T &y() requires(N > 1) { return data_[1]; }
	/// @brief Access z component (if N > 2)
	constexpr T &z() requires(N > 2) { return data_[2]; }
	/// @brief Access w component (if N > 3)
	constexpr T &w() requires(N > 3) { return data_[3]; }

	/// @brief Const access to x component
	constexpr const T &x() const requires(N > 0) { return data_[0]; }
	/// @brief Const access to y component
	constexpr const T &y() const requires(N > 1) { return data_[1]; }
	/// @brief Const access to z component
	constexpr const T &z() const requires(N > 2) { return data_[2]; }
	/// @brief Const access to w component
	constexpr const T &w() const requires(N > 3) { return data_[3]; }

	// Arithmetic operations

	/** @brief Vector addition */
	constexpr Vector operator+(const Vector &rhs) const
	{
		Vector result;

		for (std::size_t i = 0; i < N; ++i) {
			result[i] = data_[i] + rhs[i];
		}

		return result;
	}

	/** @brief Vector subtraction */
	constexpr Vector operator-(const Vector &rhs) const
	{
		Vector result;

		for (std::size_t i = 0; i < N; ++i) {
			result[i] = data_[i] - rhs[i];
		}

		return result;
	}

	/** @brief Scalar multiplication */
	constexpr Vector operator*(T scalar) const
	{
		Vector result;

		for (std::size_t i = 0; i < N; ++i) {
			result[i] = data_[i] * scalar;
		}

		return result;
	}

	/** @brief Scalar division */
	constexpr Vector operator/(T scalar) const
	{
		Vector result;

		for (std::size_t i = 0; i < N; ++i) {
			result[i] = data_[i] / scalar;
		}

		return result;
	}

	/**
	 * @brief Left scalar multiplication
	 * @param scalar Scalar value
	 * @param v Vector to multiply
	 * @return Result of scalar * v
	 */
	friend constexpr Vector operator*(T scalar, const Vector &v)
	{
		return v * scalar;
	}

	/**
	 * @brief Dot product
	 * @param rhs Other vector
	 * @return Dot product result
	 */
	constexpr T dot(const Vector &rhs) const
	{
		T result = T{};

		for (std::size_t i = 0; i < N; ++i) {
			result += data_[i] * rhs[i];
		}

		return result;
	}

	/** @brief Squared norm */
	constexpr T normSquared() const
	{
		return dot(*this);
	}

	/** @brief Euclidean norm */
	constexpr T norm() const
	{
		return std::sqrt(normSquared());
	}

	/** @brief Normalized vector */
	constexpr Vector normalized() const
	{
		T n = norm();
		return (n != T{0}) ? (*this / n) : *this;
	}

	/**
	 * @brief Cross product (3D only)
	 * @param rhs Other vector
	 * @return Cross product vector
	 */
	constexpr Vector cross(const Vector &rhs) const requires(N == 3)
	{
		return Vector{
			data_[1] *rhs[2] - data_[2] *rhs[1],
			data_[2] *rhs[0] - data_[0] *rhs[2],
			data_[0] *rhs[1] - data_[1] *rhs[0]
		};
	}

	/** @brief Component-wise multiplication */
	constexpr Vector cwiseProduct(const Vector &rhs) const
	{
		Vector result;

		for (std::size_t i = 0; i < N; ++i) {
			result[i] = data_[i] * rhs[i];
		}

		return result;
	}

	/** @brief Component-wise absolute value */
	constexpr Vector cwiseAbs() const
	{
		Vector result;

		for (std::size_t i = 0; i < N; ++i) {
			result[i] = std::abs(data_[i]);
		}

		return result;
	}

	/** @brief Converts to std::array */
	constexpr std::array<T, N> to_array() const
	{
		return data_;
	}
};

/** @brief Type alias for 2D vector */
template<typename T> using Vector2 = Vector<T, 2>;

/** @brief Type alias for 3D vector */
template<typename T> using Vector3 = Vector<T, 3>;

/** @brief Type alias for 4D vector */
template<typename T> using Vector4 = Vector<T, 4>;

} // namespace Cadence
