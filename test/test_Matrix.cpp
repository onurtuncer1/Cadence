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

TEST_CASE("Matrix4x4 Identity factory", "[Matrix]") {
    using Mat = Matrix4x4<float>;
    const auto mat = Mat::identity();

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

TEST_CASE("Matrix4x4 Zero factory", "[Matrix]") {
    using Mat = Matrix4x4<float>;
    const Mat z = Mat::zero();
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 4; ++j)
            REQUIRE(z(i, j) == Approx(0.0f));
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

TEST_CASE("Matrix4x4: Identity and Zero multiplication (both sides)", "[Matrix]") {
    using Mat = Matrix4x4<float>;

    // Build a deterministic test matrix with distinct values
    Mat M{};
    {
        float v = 1.0f;
        for (size_t r = 0; r < 4; ++r)
            for (size_t c = 0; c < 4; ++c)
                M(r, c) = v++;
    }

    const Mat I = Mat::identity();
    const Mat Z = Mat::zero();

    SECTION("Right-multiply by Identity preserves M") {
        const Mat R = M * I;
        for (size_t r = 0; r < 4; ++r)
            for (size_t c = 0; c < 4; ++c)
                REQUIRE(R(r, c) == Approx(M(r, c)));
    }

    SECTION("Left-multiply by Identity preserves M") {
        const Mat R = I * M;
        for (size_t r = 0; r < 4; ++r)
            for (size_t c = 0; c < 4; ++c)
                REQUIRE(R(r, c) == Approx(M(r, c)));
    }

    SECTION("Right-multiply by Zero yields Zero") {
        const Mat R = M * Z;
        for (size_t r = 0; r < 4; ++r)
            for (size_t c = 0; c < 4; ++c)
                REQUIRE(R(r, c) == Approx(0.0f));
    }

    SECTION("Left-multiply by Zero yields Zero") {
        const Mat R = Z * M;
        for (size_t r = 0; r < 4; ++r)
            for (size_t c = 0; c < 4; ++c)
                REQUIRE(R(r, c) == Approx(0.0f));
    }
}