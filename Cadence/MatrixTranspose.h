

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace Cadence {

template<std::floating_point T, size_t Rows, size_t Cols>
constexpr auto Transpose(const std::array<std::array<T, Cols>, Rows>& m) {
    std::array<std::array<T, Rows>, Cols> result{};
    for (size_t i = 0; i < Rows; ++i)
        for (size_t j = 0; j < Cols; ++j)
            result[j][i] = m[i][j];
    return result;
}

}