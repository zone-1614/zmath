#pragma once

#include <limits>

namespace zmath {

const double epsilon = 1e-12;

const double inf = std::numeric_limits<double>::infinity();

const double ninf = -1 * inf; // negative inf

const int zinf = std::numeric_limits<int>::max(); // z 代表 整数Z

const int nzinf = std::numeric_limits<int>::min();

const double pi = 3.14159265358979323846;

const double exp = 2.71828182845904523536;

bool is_zero(double d) {
    return std::abs(d) < epsilon;
}

bool eq(double a, double b) {
    return std::abs(a - b) < epsilon;
}

constexpr double deg2rad(double deg) {
    return deg * pi / 180;
}

constexpr double rad2deg(double rad) {
    return 180 * rad / pi;
}

}