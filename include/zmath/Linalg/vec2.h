#pragma once

#include <iostream>

#include <fmt/format.h>

#include <zmath/utils.h>

namespace zmath {

class Vec2 {
public:
    double x {0.0};
    double y {0.0};
    Vec2() = default;
    Vec2(double x, double y) : x(x), y(y) { }

    // 点乘
    double dot(const Vec2& rhs) const {
        return x * rhs.x + y * rhs.y;
    }

    // 范数
    double norm() const {
        return std::hypot(x, y);
    }

    // 距离
    double distance(const Vec2& v) const {
        return std::hypot(x - v.x, y - v.y);
    } 

    // 单位化
    Vec2 normalize() const {
        Vec2 v = *this;
        v /= v.norm();
        return v;
    }

    // 两个向量之间的夹角, 返回弧度
    double angle(const Vec2& v) const {
        double d = dot(v);
        d /= norm();
        d /= v.norm();
        // clamp
        d = std::max(std::min(d, 1.0), -1.0);
        return std::acos(d);
    }

    Vec2 project(const Vec2& v) const {
        const double p = norm() * std::cos(angle(v));
        return v.normalize() *= p;
    }

    Vec2 operator+(const Vec2& rhs) {
        Vec2 ret = *this;
        ret.x += rhs.x;
        ret.y += rhs.y;
        return ret;
    }

    Vec2 operator-(const Vec2& rhs) {
        Vec2 ret = *this;
        ret.x -= rhs.x;
        ret.y -= rhs.y;
        return ret;
    }

    Vec2& operator+=(const Vec2& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    Vec2& operator-=(const Vec2& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    Vec2& operator*=(double scaler) {
        x *= scaler;
        y *= scaler;
        return *this;
    }

    Vec2& operator/=(double scaler) {
        x /= scaler;
        y /= scaler;
        return *this;
    }

    bool operator==(const Vec2& rhs) const {
        return zmath::eq(x, rhs.x) && zmath::eq(y, rhs.y);
    }

    bool operator!=(const Vec2& rhs) const {
        return !(*this == rhs);
    }

    std::string to_string() const {
        return fmt::format("Vec2 ( {:.3f}, {:.3f} )", x, y);
    }

    friend std::ostream& operator<<(std::ostream& os, const Vec2& v) {
        os << v.to_string();
        return os;
    }

    void print() const {
        std::cout << *this << std::endl;
    }

    static Vec2 up() {
        return Vec2(0.0, 1.0);
    }

    static Vec2 down() {
        return Vec2(0.0, -1.0);
    }

    static Vec2 left() {
        return Vec2(-1.0, 0.0);
    }

    static Vec2 right() {
        return Vec2(1.0, 0.0);
    }
};

}