//------------------------------------------------------------------------------
// Project: Cadence Math Library
// Copyright (c) 2025, Onur Tuncer, PhD, Istanbul Technical University
// SPDX-License-Identifier: MIT
// License-Filename: LICENSE
//------------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "Cadence/Transform.h"

using Catch::Approx;

namespace {
template <typename T>
inline void expect_vec3_eq(const Cadence::Vector3<T>& a,
                           const Cadence::Vector3<T>& b,
                           T absTol = static_cast<T>(1e-6),
                           T relTol = static_cast<T>(1e-6)) {
    using Catch::Approx;

    // Use both absolute and relative tolerances; absolute handles near-zero.
    REQUIRE(a[0] == Approx(b[0]).margin(absTol).epsilon(relTol));
    REQUIRE(a[1] == Approx(b[1]).margin(absTol).epsilon(relTol));
    REQUIRE(a[2] == Approx(b[2]).margin(absTol).epsilon(relTol));
}
} // namespace

TEST_CASE("Transform: identity behavior", "[Transform]") {
    using T  = float;
    using TF = Cadence::Transform<T>;
    using V3 = Cadence::Vector3<T>;

    const TF I{};  // default ctor -> identity
    const V3 p{1, 2, 3};
    const V3 v{4, 5, 6};

    SECTION("identity leaves point unchanged") {
        expect_vec3_eq(I.transform_point(p), p);
    }

    SECTION("identity leaves vector unchanged") {
        expect_vec3_eq(I.transform_vector(v), v);
    }
}

TEST_CASE("Transform: translation", "[Transform]") {
    using T  = float;
    using TF = Cadence::Transform<T>;
    using V3 = Cadence::Vector3<T>;

    const TF Tr = TF::translation(10, -2, 0.5f);

    SECTION("translation affects points") {
        const V3 p{1, 2, 3};
        expect_vec3_eq(Tr.transform_point(p), V3{11, 0, 3.5f});
    }

    SECTION("translation does not affect vectors") {
        const V3 v{4, 5, 6};
        expect_vec3_eq(Tr.transform_vector(v), v);
    }
}

TEST_CASE("Transform: scaling", "[Transform]") {
    using T  = float;
    using TF = Cadence::Transform<T>;
    using V3 = Cadence::Vector3<T>;

    const TF S = TF::scaling(2, 3, -4);

    SECTION("scaling affects points") {
        const V3 p{1, -2, 0.5f};
        expect_vec3_eq(S.transform_point(p), V3{2, -6, -2});
    }

    SECTION("scaling affects vectors") {
        const V3 v{1, 2, 3};
        expect_vec3_eq(S.transform_vector(v), V3{2, 6, -12});
    }
}

TEST_CASE("Transform: axis-aligned rotations (π/2)", "[Transform]") {
    using T  = float;
    using TF = Cadence::Transform<T>;
    using V3 = Cadence::Vector3<T>;

    constexpr T PI = static_cast<T>(3.14159265358979323846);
    constexpr T H  = PI / static_cast<T>(2); // 90 degrees

    SECTION("rotation_x") {
        const TF Rx = TF::rotation_x(H);
        expect_vec3_eq(Rx.transform_vector(V3{0, 1, 0}), V3{0, 0, 1});
        expect_vec3_eq(Rx.transform_vector(V3{0, 0, 1}), V3{0,-1, 0});
    }

    SECTION("rotation_y") {
        const TF Ry = TF::rotation_y(H);
        expect_vec3_eq(Ry.transform_vector(V3{0, 0, 1}), V3{1, 0, 0});
        expect_vec3_eq(Ry.transform_vector(V3{1, 0, 0}), V3{0, 0,-1});
    }

    SECTION("rotation_z") {
        const TF Rz = TF::rotation_z(H);
        expect_vec3_eq(Rz.transform_vector(V3{1, 0, 0}), V3{0, 1, 0});
        expect_vec3_eq(Rz.transform_vector(V3{0, 1, 0}), V3{-1,0, 0});
    }
}

TEST_CASE("Transform: rotation around arbitrary axis", "[Transform]") {
    using T  = float;
    using TF = Cadence::Transform<T>;
    using V3 = Cadence::Vector3<T>;

    constexpr T PI = static_cast<T>(3.14159265358979323846);
    constexpr T H  = PI / static_cast<T>(2); // 90 degrees

    // Use non-unit axis to verify internal normalization
    const TF Rz_from_axis = TF::rotation_axis(V3{0, 0, 5}, H);

    SECTION("equivalence with rotation_z for 90° about Z") {
        expect_vec3_eq(Rz_from_axis.transform_vector(V3{1, 0, 0}), V3{0, 1, 0});
        expect_vec3_eq(Rz_from_axis.transform_vector(V3{0, 1, 0}), V3{-1, 0, 0});
    }

    SECTION("vector parallel to axis remains on axis") {
        // Any rotation around Z keeps Z-component direction
        expect_vec3_eq(Rz_from_axis.transform_vector(V3{0, 0, 3}), V3{0, 0, 3});
    }
}

TEST_CASE("Transform: composition order", "[Transform]") {
    using T  = float;
    using TF = Cadence::Transform<T>;
    using V3 = Cadence::Vector3<T>;

    constexpr T PI = static_cast<T>(3.14159265358979323846);
    constexpr T H  = PI / static_cast<T>(2); // 90 degrees

    const TF Rz = TF::rotation_z(H);     // rotate +90° about Z
    const TF Tx = TF::translation(1, 0, 0);

    // Composition is left-to-right (result.m = this->m * other.m):
    // (Rz * Tx)(p) = Rz( Tx(p) )
    SECTION("Rz * Tx then applied to origin -> (0,1,0)") {
        const TF C = Rz * Tx;
        expect_vec3_eq(C.transform_point(V3{0,0,0}), V3{0,1,0});
    }

    SECTION("Tx * Rz then applied to origin -> (1,0,0)") {
        const TF C = Tx * Rz;
        expect_vec3_eq(C.transform_point(V3{0,0,0}), V3{1,0,0});
    }
}

TEST_CASE("Transform: identity is neutral in composition", "[Transform]") {
    using T  = float;
    using TF = Cadence::Transform<T>;
    using V3 = Cadence::Vector3<T>;

    const TF I{};
    const TF T1 = TF::translation(3, -2, 5);

    SECTION("T1 * I == T1") {
        const TF C = T1 * I;
        expect_vec3_eq(C.transform_point(V3{1,2,3}), T1.transform_point(V3{1,2,3}));
    }

    SECTION("I * T1 == T1") {
        const TF C = I * T1;
        expect_vec3_eq(C.transform_point(V3{1,2,3}), T1.transform_point(V3{1,2,3}));
    }
}
