//------------------------------------------------------------------------------
// Project: Cadence Math Library
// Copyright (c) 2025, Onur Tuncer, PhD, Istanbul Technical University
// SPDX-License-Identifier: MIT
// License-Filename: LICENSE
//------------------------------------------------------------------------------

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <concepts>
#include <type_traits>

#include "Cadence/Matrix.h"

/**
 * @file Transform.h
 * @brief Defines the Cadence::Transform class for affine 3D transformations.
 */

namespace Cadence
{

/**
 * @class Transform
 * @brief Represents a 4x4 affine transformation matrix in 3D space.
 *
 * Provides factory functions for translation, scaling, rotation around axes,
 * and rotation around an arbitrary axis using the Rodrigues formula.
 *
 * @tparam T Floating-point type used for internal representation.
 */
template<typename T>
class Transform
{
	Matrix4x4<T> m{}; ///< Internal 4x4 transformation matrix

public:
	/**
	 * @brief Default constructor initializes to identity transformation.
	 */
	constexpr Transform()
	{
		for (size_t i = 0; i < 4; ++i) {
			m(i, i) = 1;
		}
	}

	/**
	 * @brief Constructs a translation transformation.
	 * @param x X-axis translation
	 * @param y Y-axis translation
	 * @param z Z-axis translation
	 * @return A translation transform
	 */
	constexpr static Transform translation(T x, T y, T z)
	{
		Transform t;
		t.m(0, 3) = x;
		t.m(1, 3) = y;
		t.m(2, 3) = z;
		return t;
	}

	/**
	 * @brief Constructs a scaling transformation.
	 * @param x Scale along X
	 * @param y Scale along Y
	 * @param z Scale along Z
	 * @return A scaling transform
	 */
	constexpr static Transform scaling(T x, T y, T z)
	{
		Transform t;
		t.m(0, 0) = x;
		t.m(1, 1) = y;
		t.m(2, 2) = z;
		return t;
	}

	/**
	 * @brief Constructs a rotation about the X-axis.
	 * @param angle_rad Angle in radians
	 * @return A transform representing X-axis rotation
	 */
	constexpr static Transform rotation_x(T angle_rad)
	{
		Transform t;
		T c = std::cos(angle_rad);
		T s = std::sin(angle_rad);
		t.m(1, 1) = c;
		t.m(1, 2) = -s;
		t.m(2, 1) = s;
		t.m(2, 2) = c;
		return t;
	}

	/**
	 * @brief Constructs a rotation about the Y-axis.
	 * @param angle_rad Angle in radians
	 * @return A transform representing Y-axis rotation
	 */
	constexpr static Transform rotation_y(T angle_rad)
	{
		Transform t;
		T c = std::cos(angle_rad);
		T s = std::sin(angle_rad);
		t.m(0, 0) = c;
		t.m(0, 2) = s;
		t.m(2, 0) = -s;
		t.m(2, 2) = c;
		return t;
	}

	/**
	 * @brief Constructs a rotation about the Z-axis.
	 * @param angle_rad Angle in radians
	 * @return A transform representing Z-axis rotation
	 */
	constexpr static Transform rotation_z(T angle_rad)
	{
		Transform t;
		T c = std::cos(angle_rad);
		T s = std::sin(angle_rad);
		t.m(0, 0) = c;
		t.m(0, 1) = -s;
		t.m(1, 0) = s;
		t.m(1, 1) = c;
		return t;
	}

	/**
	 * @brief Constructs a rotation around an arbitrary axis.
	 * @param axis Rotation axis (must be non-zero)
	 * @param angle_rad Rotation angle in radians
	 * @return A transform representing the rotation
	 */
	constexpr static Transform rotation_axis(Vector3<T> axis, T angle_rad)
	{
		axis = normalize(axis);
		T c = std::cos(angle_rad);
		T s = std::sin(angle_rad);
		T t = 1 - c;

		Transform rot;
		rot.m(0, 0) = t * axis[0] * axis[0] + c;
		rot.m(0, 1) = t * axis[0] * axis[1] - s * axis[2];
		rot.m(0, 2) = t * axis[0] * axis[2] + s * axis[1];

		rot.m(1, 0) = t * axis[0] * axis[1] + s * axis[2];
		rot.m(1, 1) = t * axis[1] * axis[1] + c;
		rot.m(1, 2) = t * axis[1] * axis[2] - s * axis[0];

		rot.m(2, 0) = t * axis[0] * axis[2] - s * axis[1];
		rot.m(2, 1) = t * axis[1] * axis[2] + s * axis[0];
		rot.m(2, 2) = t * axis[2] * axis[2] + c;

		return rot;
	}

	/**
	 * @brief Composes two transformations via matrix multiplication.
	 * @param other Right-hand transform
	 * @return Resulting composed transform
	 */
	constexpr Transform operator*(const Transform &other) const
	{
		Transform result;

		for (size_t i = 0; i < 4; ++i) {
			for (size_t j = 0; j < 4; ++j) {
				T sum = 0;

				for (size_t k = 0; k < 4; ++k) {
					sum += m(i, k) * other.m(k, j);
				}

				result.m(i, j) = sum;
			}
		}

		return result;
	}

	/**
	 * @brief Transforms a 3D point (implicitly sets w=1).
	 * @param p 3D point to transform
	 * @return Transformed point
	 */
	constexpr Vector3<T> transform_point(const Vector3<T> &p) const
	{
		Vector4<T> homog{p[0], p[1], p[2], 1};
		homog = m * homog;
		return {homog[0], homog[1], homog[2]};
	}

	/**
	 * @brief Transforms a 3D vector (implicitly sets w=0).
	 * @param v 3D direction vector
	 * @return Transformed vector
	 */
	constexpr Vector3<T> transform_vector(const Vector3<T> &v) const
	{
		Vector4<T> homog{v[0], v[1], v[2], 0};
		homog = m * homog;
		return {homog[0], homog[1], homog[2]};
	}

private:
	/**
	 * @brief Normalizes a 3D vector.
	 * @param v Input vector
	 * @return Normalized vector
	 */
	constexpr static Vector3<T> normalize(Vector3<T> v)
	{
		T len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		return {v[0] / len, v[1] / len, v[2] / len};
	}
};

} // namespace Cadence