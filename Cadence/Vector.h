//============================================================================//
// Cadence - Constexpr Math for CAD Kernels
// SPDX-License-Identifier: MIT
// Author: Onur Tuncer, PhD onur.tuncer@itu.edu.tr
// File: [filename]
// Copyright (c) 2025 Onur Tuncer, PhD
// Desc: [1-2 line description of purpose]
//============================================================================//

#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>

namespace Cadence {

template<std::floating_point T, std::size_t N>
class Vector {
    std::array<T, N> data_{};

public:
    static constexpr std::size_t static_size = N;

    // Constructors
    constexpr Vector() = default;

    template<typename... Args>
        requires (sizeof...(Args) == N) && (std::conjunction_v<std::is_convertible<Args, T>...>)
    constexpr Vector(Args... args) : data_{static_cast<T>(args)...} {}

    constexpr Vector(std::array<T, N> init) : data_{init} {}

    // Accessors
    constexpr T& operator[](std::size_t i) { return data_[i]; }
    constexpr const T& operator[](std::size_t i) const { return data_[i]; }

    constexpr std::size_t size() const { return N; }

    // Named component access (like Eigen)
    constexpr T& x() requires (N > 0) { return data_[0]; }
    constexpr T& y() requires (N > 1) { return data_[1]; }
    constexpr T& z() requires (N > 2) { return data_[2]; }
    constexpr T& w() requires (N > 3) { return data_[3]; }

    constexpr const T& x() const requires (N > 0) { return data_[0]; }
    constexpr const T& y() const requires (N > 1) { return data_[1]; }
    constexpr const T& z() const requires (N > 2) { return data_[2]; }
    constexpr const T& w() const requires (N > 3) { return data_[3]; }

    // Arithmetic
    constexpr Vector operator+(const Vector& rhs) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = data_[i] + rhs[i];
        return result;
    }

    constexpr Vector operator-(const Vector& rhs) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = data_[i] - rhs[i];
        return result;
    }

    constexpr Vector operator*(T scalar) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = data_[i] * scalar;
        return result;
    }

    constexpr Vector operator/(T scalar) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = data_[i] / scalar;
        return result;
    }

    // Scalar multiplication from left
    friend constexpr Vector operator*(T scalar, const Vector& v) {
        return v * scalar;
    }

    // Dot product
    constexpr T dot(const Vector& rhs) const {
        T result = T{};
        for (std::size_t i = 0; i < N; ++i)
            result += data_[i] * rhs[i];
        return result;
    }

    // Norm and normalization
    constexpr T normSquared() const {
        return dot(*this);
    }

    constexpr T norm() const {
        return std::sqrt(normSquared());
    }

    constexpr Vector normalized() const {
        T n = norm();
        return (n != T{0}) ? (*this / n) : *this;
    }

    // Cross product (3D only)
    constexpr Vector cross(const Vector& rhs) const requires (N == 3) {
        return Vector{
            data_[1] * rhs[2] - data_[2] * rhs[1],
            data_[2] * rhs[0] - data_[0] * rhs[2],
            data_[0] * rhs[1] - data_[1] * rhs[0]
        };
    }

    // Component-wise operations
    constexpr Vector cwiseProduct(const Vector& rhs) const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = data_[i] * rhs[i];
        return result;
    }

    constexpr Vector cwiseAbs() const {
        Vector result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = std::abs(data_[i]);
        return result;
    }

    constexpr std::array<T, N> to_array() const {
        return data_;
    }
};

// Aliases
template<typename T> using Vector2 = Vector<T, 2>;
template<typename T> using Vector3 = Vector<T, 3>;
template<typename T> using Vector4 = Vector<T, 4>;

} // namespace Cadence
