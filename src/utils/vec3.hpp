#ifndef ANALYSISGROUPC_UTILS_VEC3_HPP_
#define ANALYSISGROUPC_UTILS_VEC3_HPP_

#include <cmath>
#include <iostream>

#include <TVector3.h>

struct vec3 {

    double x, y, z;

    constexpr vec3() noexcept : x(0.0), y(0.0), z(0.0) {}
    constexpr vec3(double x_, double y_, double z_) noexcept : x(x_), y(y_), z(z_) {}

    explicit vec3(const TVector3& v) : x(v.X()), y(v.Y()), z(v.Z()) {}
    
    constexpr vec3(const vec3& v) noexcept = default;
    constexpr vec3(vec3&& v) noexcept = default;

    ~vec3() noexcept = default;

    constexpr vec3& operator=(const vec3& v) noexcept = default;
    constexpr vec3& operator=(vec3&& v) noexcept = default;
    
    vec3& operator=(const TVector3& v) {
        x = v.X();
        y = v.Y();
        z = v.Z();
        return *this;
    }

    constexpr vec3& operator+=(const vec3& v) noexcept {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    constexpr vec3& operator-=(const vec3& v) noexcept {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    constexpr vec3 operator-() const noexcept {
        return vec3{-x, -y, -z};
    }

    constexpr vec3& operator*=(double s) noexcept {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    constexpr vec3& operator/=(double s) noexcept {
        return *this *= (1.0 / s);
    }

};

inline constexpr vec3 operator+(const vec3& lhs, const vec3& rhs) noexcept {
    return vec3{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

inline constexpr vec3 operator-(const vec3& lhs, const vec3& rhs) noexcept {
    return vec3{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}


inline constexpr vec3 operator*(const vec3& v, double s) noexcept {
    return vec3{v.x * s, v.y * s, v.z * s};
}

inline constexpr vec3 operator*(double s, const vec3& v) noexcept {
    return v * s;
}

inline constexpr vec3 operator/(const vec3& v, double s) noexcept {
    return v * (1.0 / s);
}

inline constexpr vec3 cross(const vec3& lhs, const vec3& rhs) noexcept {
    return vec3{lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x};
}

inline constexpr double dot(const vec3& lhs, const vec3& rhs) noexcept {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
} 

inline constexpr double mag2(const vec3& v) noexcept {
    return dot(v, v);
}

inline constexpr double mag(const vec3& v) noexcept {
    return std::sqrt(mag2(v));
}

inline constexpr vec3 unit(const vec3& v) noexcept {
    return v / mag(v);
}

inline constexpr double angle(const vec3& lhs, const vec3& rhs) noexcept {
    return std::acos(dot(lhs, rhs) / std::sqrt(mag2(lhs) * mag2(rhs)));
}

inline constexpr double theta(const vec3& v) noexcept {
    return std::acos(v.z / mag(v));
}

inline constexpr double phi(const vec3& v) noexcept {
    return std::atan2(v.y, v.x);
}

inline constexpr vec3 from_spherical(double r, double theta, double phi) noexcept {
    r = std::abs(r);
    double x = r * std::sin(theta) * std::cos(phi);
    double y = r * std::sin(theta) * std::sin(phi);
    double z = r * std::cos(theta);
    return vec3{x, y, z};
}

template<class _Char, class _Traits>
std::basic_ostream<_Char, _Traits>& operator<<(std::basic_ostream<_Char, _Traits>& os, const vec3& vec) {
    return os << '(' << vec.x << ", " << vec.y << ", " << vec.z << ')';
}

#endif // ANALYSISGROUPC_UTILS_VEC3_HPP_