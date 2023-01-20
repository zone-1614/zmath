#pragma once

#include <iostream>
#include <vector>

#include <zmath/fft.h>

namespace zmath {

class Polynomial {
public:
    /**
     * @brief 构造一个多项式
     * 
     * @param coef 多项式的系数. 例如: 多项式是  3x^2 - 2x + 1, 则 coef 是 { 3, -2, 1 } 
     */
    Polynomial(const std::vector<double>& coef) : coef_(coef) {

    }

    /**
     * @brief 返回多项式的度
     */
    size_t deg() const {
        return coef_.size() - 1;
    }

    /**
     * @brief 返回多项式的系数
     */
    std::vector<double> coef() const {
        return coef_;
    }

    /**
     * @brief 返回多项式的导数
     */
    Polynomial derivative() const {
        if (coef_.size() <= 1) return Polynomial({0});
    }

    /**
     * @brief 首一  TODO
     */
    Polynomial monic() const {
        return Polynomial({0});
    }

    // TODO
    double operator()(double x) const {
        return 0.0f;
    }

    // TODO
    Polynomial operator+(const Polynomial& rhs) const {
        return Polynomial({0});
    }

    // TODO
    Polynomial operator-(const Polynomial& rhs) const {
        return Polynomial({0});
    }

    // TODO
    Polynomial operator*(const Polynomial& rhs) const {
        return Polynomial({0});
    }

    // TODO
    Polynomial operator^(const Polynomial& rhs) const {
        return Polynomial({0});
    }

    // TODO
    Polynomial& operator+=(const Polynomial& rhs) {
        return *this;
    }

    // TODO
    Polynomial& operator-=(const Polynomial& rhs) {
        return *this;
    }

    // TODO
    Polynomial& operator*=(const Polynomial& rhs) {
        return *this;
    }

    // TODO
    Polynomial& operator^=(const Polynomial& rhs) {
        return *this;
    }

    // TODO
    friend std::ostream& operator<<(std::ostream& os, const Polynomial& poly) {
        
        return os;
    }

    // TODO
    void print() const {
        std::cout << *this << std::endl;
    }

private:
    std::vector<double> coef_;
};

}