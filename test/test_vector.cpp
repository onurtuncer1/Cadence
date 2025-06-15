

#include <catch2/catch_test_macros.hpp>
#include "Cadence/Vector.h"

using Cadence::Vector3;

TEST_CASE("Vector addition works", "[vector]") {
    constexpr Vector3<double> a{1.0, 2.0, 3.0};
    constexpr Vector3<double> b{4.0, 5.0, 6.0};
    constexpr auto c = a + b;
    REQUIRE(c.x() == 5.0);
    REQUIRE(c.y() == 7.0);
    REQUIRE(c.z() == 9.0);
}
