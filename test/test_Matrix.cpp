#include <catch2/catch_test_macros.hpp>
#include "Cadence/Matrix.h"
#include "Cadence/Vector.h"

using namespace Cadence;

TEST_CASE("Matrix element access and construction", "[matrix]") {
    constexpr Matrix2x2<double> m = Matrix2x2<double>() << 1.0, 2.0, 3.0, 4.0;
    STATIC_REQUIRE(m(0, 0) == 1.0);
    STATIC_REQUIRE(m(0, 1) == 2.0);
    STATIC_REQUIRE(m(1, 0) == 3.0);
    STATIC_REQUIRE(m(1, 1) == 4.0);
}

TEST_CASE("Matrix identity generation", "[matrix]") {
    constexpr auto id = Matrix3x3<double>::identity();
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            if (i == j)
                STATIC_REQUIRE(id(i, j) == 1.0);
            else
                STATIC_REQUIRE(id(i, j) == 0.0);
        }
    }
}

TEST_CASE("Matrix * Vector multiplication", "[matrix][vector]") {
    constexpr Matrix3x3<double> m = Matrix3x3<double>() << 
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0;

    constexpr Vector3<double> v{1.0, 0.0, -1.0};
    constexpr auto result = m * v;

    STATIC_REQUIRE(result[0] == Approx(-2.0));
    STATIC_REQUIRE(result[1] == Approx(-2.0));
    STATIC_REQUIRE(result[2] == Approx(-2.0));
}

TEST_CASE("Matrix * Matrix multiplication", "[matrix]") {
    constexpr Matrix2x2<double> a = Matrix2x2<double>() << 1.0, 2.0, 3.0, 4.0;
    constexpr Matrix2x2<double> b = Matrix2x2<double>() << 5.0, 6.0, 7.0, 8.0;

    constexpr auto c = a * b;

    STATIC_REQUIRE(c(0, 0) == Approx(19.0));
    STATIC_REQUIRE(c(0, 1) == Approx(22.0));
    STATIC_REQUIRE(c(1, 0) == Approx(43.0));
    STATIC_REQUIRE(c(1, 1) == Approx(50.0));
}

TEST_CASE("Matrix scalar operations", "[matrix]") {
    constexpr Matrix2x2<double> m = Matrix2x2<double>() << 1.0, 2.0, 3.0, 4.0;

    constexpr auto scaled = m * 2.0;
    STATIC_REQUIRE(scaled(0, 0) == 2.0);
    STATIC_REQUIRE(scaled(1, 1) == 8.0);

    constexpr auto divided = m / 2.0;
    STATIC_REQUIRE(divided(0, 0) == 0.5);
    STATIC_REQUIRE(divided(1, 1) == 2.0);
}

