#pragma once

#include <iostream>

#include <fmt/format.h>

#include <zmath/utils.h>

namespace zmath {

class Vec3 {
public:
    double x {0.0};
    double y {0.0};
    double z {0.0};
    Vec3() = default;
    Vec3(double x, double y, double z) : x(x), y(y), z(z) { }

    // 点乘
    double dot(const Vec3& rhs) const {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    // 范数
    double norm() const {
        return std::hypot(x, y, z);
    }

    // 距离
    double distance(const Vec3& v) const {
        return std::hypot(x - v.x, y - v.y, z - v.z);
    } 

    // 单位化
    Vec3 normalize() const {
        Vec3 v = *this;
        v /= v.norm();
        return v;
    }

    // 两个向量之间的夹角, 返回弧度
    double angle(const Vec3& v) const {
        double d = dot(v);
        d /= norm();
        d /= v.norm();
        // clamp
        d = std::max(std::min(d, 1.0), -1.0);
        return std::acos(d);
    }

    Vec3 project(const Vec3& v) const {
        const double p = norm() * std::cos(angle(v));
        return v.normalize() *= p;
    }

    Vec3 operator+(const Vec3& rhs) {
        Vec3 ret = *this;
        ret.x += rhs.x;
        ret.y += rhs.y;
        ret.z += rhs.z;
        return ret;
    }

    Vec3 operator-(const Vec3& rhs) {
        Vec3 ret = *this;
        ret.x -= rhs.x;
        ret.y -= rhs.y;
        ret.z -= rhs.z;
        return ret;
    }

    Vec3& operator+=(const Vec3& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    Vec3& operator*=(double scaler) {
        x *= scaler;
        y *= scaler;
        z *= scaler;
        return *this;
    }

    Vec3& operator/=(double scaler) {
        x /= scaler;
        y /= scaler;
        z /= scaler;
        return *this;
    }

    bool operator==(const Vec3& rhs) const {
        return zmath::eq(x, rhs.x) && zmath::eq(y, rhs.y) && zmath::eq(z, rhs.z);
    }

    bool operator!=(const Vec3& rhs) const {
        return !(*this == rhs);
    }

    std::string to_string() const {
        return fmt::format("Vec3 ( {:.3f}, {:.3f}, {:.3f} )", x, y, z);
    }

    friend std::ostream& operator<<(std::ostream& os, const Vec3& v) {
        os << v.to_string();
        return os;
    }

    void print() const {
        std::cout << *this << std::endl;
    }

    static Vec3 forward() {
        return Vec3(0.0, 0.0, 1.0);
    }

    static Vec3 back() {
        return Vec3(0.0, 0.0, -1.0);
    }

    static Vec3 up() {
        return Vec3(0.0, 1.0, 0.0);
    } 

    static Vec3 down() {
        return Vec3(0.0, -1.0, 0.0);
    }

    static Vec3 left() {
        return Vec3(-1.0, 0.0, 0.0);
    }

    static Vec3 right() {
        return Vec3(1.0, 0.0, 0.0);
    }
    
};

}