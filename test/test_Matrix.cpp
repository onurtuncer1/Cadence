//------------------------------------------------------------------------------
// Project: Cadence Math Library
// Copyright (c) 2025, Onur Tuncer, PhD, Istanbul Technical University
// SPDX-License-Identifier: MIT
// License-Filename: LICENSE
//------------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "Cadence/Matrix.h"  // adjust path as needed

using Cadence::Matrix4x4;
using Catch::Approx;

TEST_CASE("Matrix4x4 Default Constructor", "[Matrix]")
{
	Matrix4x4<float> mat;

	for (size_t row = 0; row < 4; ++row) {
		for (size_t col = 0; col < 4; ++col) {
			if (row == col) {
				REQUIRE(mat(row, col) == Approx(1.0f));

			} else {
				REQUIRE(mat(row, col) == Approx(0.0f));
			}
		}
	}
}

TEST_CASE("Matrix4x4 Assignment and Access", "[Matrix]")
{
	Matrix4x4<float> mat;
	mat(1, 2) = 5.5f;

	REQUIRE(mat(1, 2) == Approx(5.5f));
}

TEST_CASE("Matrix4x4 Copy Constructor", "[Matrix]")
{
	Matrix4x4<float> mat1;
	mat1(2, 3) = 7.0f;

	Matrix4x4<float> mat2 = mat1;
	REQUIRE(mat2(2, 3) == Approx(7.0f));
}

TEST_CASE("Matrix4x4 Multiplication by Identity", "[Matrix]")
{
	Matrix4x4<float> identity;
	Matrix4x4<float> mat;

	// Fill with some values
	float val = 1.0f;

	for (size_t r = 0; r < 4; ++r)
		for (size_t c = 0; c < 4; ++c) {
			mat(r, c) = val++;
		}

	Matrix4x4<float> result = mat * identity;

	for (size_t r = 0; r < 4; ++r)
		for (size_t c = 0; c < 4; ++c) {
			REQUIRE(result(r, c) == Approx(mat(r, c)));
		}
}

