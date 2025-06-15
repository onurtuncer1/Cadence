#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <concepts>
#include <type_traits>

#include "Cadence/Matrix.h"

namespace Cadence{

template<typename T>
class Transform {
    Matrix4x4<T> m{};
    
public:
    constexpr Transform() {
        // Identity matrix
        for (size_t i = 0; i < 4; ++i) {
            m(i, i) = 1;
        }
    }
    
    // Translation
    constexpr static Transform translation(T x, T y, T z) {
        Transform t;
        t.m(0, 3) = x;
        t.m(1, 3) = y;
        t.m(2, 3) = z;
        return t;
    }
    
    // Scaling
    constexpr static Transform scaling(T x, T y, T z) {
        Transform t;
        t.m(0, 0) = x;
        t.m(1, 1) = y;
        t.m(2, 2) = z;
        return t;
    }
    
    // Rotation around X-axis
    constexpr static Transform rotation_x(T angle_rad) {
        Transform t;
        T c = std::cos(angle_rad);
        T s = std::sin(angle_rad);
        t.m(1, 1) = c;
        t.m(1, 2) = -s;
        t.m(2, 1) = s;
        t.m(2, 2) = c;
        return t;
    }
    
    // Rotation around Y-axis
    constexpr static Transform rotation_y(T angle_rad) {
        Transform t;
        T c = std::cos(angle_rad);
        T s = std::sin(angle_rad);
        t.m(0, 0) = c;
        t.m(0, 2) = s;
        t.m(2, 0) = -s;
        t.m(2, 2) = c;
        return t;
    }
    
    // Rotation around Z-axis
    constexpr static Transform rotation_z(T angle_rad) {
        Transform t;
        T c = std::cos(angle_rad);
        T s = std::sin(angle_rad);
        t.m(0, 0) = c;
        t.m(0, 1) = -s;
        t.m(1, 0) = s;
        t.m(1, 1) = c;
        return t;
    }
    
    // Arbitrary axis rotation (Rodrigues formula)
    constexpr static Transform rotation_axis(Vector3<T> axis, T angle_rad) {
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
    
    // Matrix multiplication (transform composition)
    constexpr Transform operator*(const Transform& other) const {
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
    
    // Transform a 3D point (adds implicit w=1)
    constexpr Vector3<T> transform_point(const Vector3<T>& p) const {
        Vector4<T> homog{p[0], p[1], p[2], 1};
        homog = m * homog;
        return {homog[0], homog[1], homog[2]};
    }
    
    // Transform a 3D vector (adds implicit w=0)
    constexpr Vector3<T> transform_vector(const Vector3<T>& v) const {
        Vector4<T> homog{v[0], v[1], v[2], 0};
        homog = m * homog;
        return {homog[0], homog[1], homog[2]};
    }
    
private:
    constexpr static Vector3<T> normalize(Vector3<T> v) {
        T len = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        return {v[0]/len, v[1]/len, v[2]/len};
    }
};

}