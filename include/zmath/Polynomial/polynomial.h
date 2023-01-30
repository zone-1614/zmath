#pragma once

#include <iostream>
#include <vector>
#include <float.h>
#include <algorithm>

#include <zmath/utils/fft.h>

namespace zmath {

class Polynomial {
public:
    Polynomial() {
        coef_.push_back(0);
    }

    Polynomial(size_t deg) {
        for (int i = 0; i <= deg; i++) {
            coef_.push_back(0);
        }
    }

    /**
     * @brief 构造一个多项式
     * 
     * @param coef 多项式的系数. 例如: 多项式是  3x^2 - 2x + 1, 则 coef 是 { 3, -2, 1 } 
     */
    Polynomial(const std::vector<double>& coef) {
        // 去掉前面的0, 比如传入的参数可能为 { 0, 0, 1, 2, 3 }
        auto first_not_zero = std::find_if_not(coef.begin(), coef.end(), [](auto d) {
            return is_zero(d);
        });
        if (first_not_zero == coef.end()) {
            // 输入参数全是0
            coef_.push_back(0);
        } else {
            for (auto it = first_not_zero; it != coef.end(); it++) {
                coef_.push_back(*it);
            }
        }
    }

    /**
     * @brief 设置多项式的系数
     * 
     * @param deg 对应项的次数
     * @param val 对应项的值
     */
    void set_coef(int deg, double val) {
        if (deg < 0) return ;
        if (deg > this->deg()) return ;
        int p = this->deg() - deg; // position
        if (is_zero(val)) {
            coef_[p] = 0;
        } else {
            coef_[p] = val;
        }
        
        if (!is_zero(coef_[0]))
            return ;

        // 防止设置第一个系数为0
        auto first_not_zero = std::find_if_not(coef_.begin(), coef_.end(), [](auto d){
            return is_zero(d);
        });
        if (first_not_zero == coef_.end()) {
            coef_.clear();
            coef_.push_back(0);
        } else {
            coef_.erase(coef_.begin(), first_not_zero);
        }
    }

    /**
     * @brief 返回多项式的度
     */
    int deg() const {
        if (coef_.size() == 0) {
            return zmath::nzinf;
        } else {
            if (coef_.size() == 1 && is_zero(coef_[0]))
                return zmath::nzinf;
            return coef_.size() - 1;
        }
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
        auto deriv = Polynomial(*this);
        const int n = deriv.deg();
        if (n < 1) return Polynomial({0}); // 只有常数项

        for (int i = 0; i < n; i++) {
            deriv.coef_[i] *= (n - i);
        }
        deriv.coef_.pop_back();
        
        return deriv;
    }

    // TODO
    /**
     * @brief 积分
     */
    Polynomial integral() const ;

    /**
     * @brief 返回首一多项式
     */
    Polynomial monic() const {
        if (deg() == zmath::nzinf) {
            return Polynomial({0});
        }

        auto monic = Polynomial(*this);
        const double d = coef_[0];

        for (auto& co : monic.coef_) {
            co /= d;
        }

        return monic;
    }

    double operator()(double x) const {
        double ret = 0.0;
        const int n = coef_.size();
        double pow = 1.0;
        for (int i = n - 1; i >= 0; i--) {
            ret += pow * coef_[i];
            pow *= x;
        }
        return ret;
    }

    Polynomial operator+(const Polynomial& rhs) const {
        int deg = std::max(this->deg(), rhs.deg());
        std::vector<double> ret_coef; // 逆序的相加结果系数
        auto it1 = this->coef_.rbegin();
        auto it2 = rhs.coef_.rbegin();

        while (it1 != this->coef_.rend() && it2 != rhs.coef_.rend()) {
            ret_coef.push_back(*it1 + *it2);
            it1++; it2++;
        }
        while (it1 != this->coef_.rend()) {
            ret_coef.push_back(*it1);
            it1++;
        }
        while (it2 != rhs.coef_.rend()) {
            ret_coef.push_back(*it2);
            it2++;
        }

        std::reverse(ret_coef.begin(), ret_coef.end());
        return Polynomial(ret_coef);
    }

    Polynomial operator-(const Polynomial& rhs) const {
        auto rhs_copy = rhs;
        for (auto& co : rhs_copy.coef_) {
            co = -co;
        }
        return *this + rhs_copy;
    }

    friend Polynomial operator*(const double& k, const Polynomial& rhs) {
        Polynomial ret = rhs;
        for (auto& co : ret.coef_) {
            co *= k;
        }
        return ret;
    }

    Polynomial operator*(const double& k) const {
        Polynomial ret = *this;
        for (auto& co : ret.coef_) {
            co *= k;
        }
        return ret;
    }

    Polynomial operator*(const Polynomial& rhs) const {
        // 计算 n=2^t 使得 n >= new_deg  刚好是2^t是因为fft.h里的fft只能运行于2^t次单位根群上
        int new_deg = deg() + rhs.deg() + 1;
        int n = 1;
        while (n < new_deg) {
            n *= 2;
        }

        // 将系数反转, 同时增长到n
        std::vector<double> coef1 = this->coef_;
        std::vector<double> coef2 = rhs.coef_;
        std::reverse(coef1.begin(), coef1.end()); 
        std::reverse(coef2.begin(), coef2.end());
        while (coef1.size() < n)
            coef1.push_back(0);
        while (coef2.size() < n)
            coef2.push_back(0);

        // 将系数转为 complex array, 进行 fft
        complex_array c1(n), c2(n);
        for (int i = 0; i < n; i++) {
            c1[i] = complex(coef1[i]);
            c2[i] = complex(coef2[i]);
        }
        fft(c1); fft(c2);
        for (int i = 0; i < n; i++) {
            c1[i] = c1[i] * c2[i];
        }
        ifft(c1);

        // 转化为结果多项式的系数的vector
        std::vector<double> ret_coef(c1.size());
        for (int i = 0; i < n; i++) {
            ret_coef[i] = c1[n - 1 - i].real();
        }

        return Polynomial(ret_coef);
    }

    Polynomial operator^(size_t t) const {
        Polynomial ret, a = *this;
        ret.coef_[0] = 1;
        if (is_zero(ret.coef_[0]))
            return ret;

        while (t) {
            if (t & 1)
                ret *= a;
            a *= a;
            t >>= 1;
        }
        
        return ret;
    }

    Polynomial& operator+=(const Polynomial& rhs) {
        *this = *this + rhs;
        return *this;
    }

    Polynomial& operator-=(const Polynomial& rhs) {
        *this = *this - rhs;
        return *this;
    }

    Polynomial& operator*=(const double& k) {
        *this = *this * k;
        return *this;
    }

    Polynomial& operator*=(const Polynomial& rhs) {
        *this = *this * rhs;
        return *this;
    }

    Polynomial& operator^=(double t) {
        *this = *this^t;
        return *this;
    }

    std::string to_string() const {
        const int n = deg();

        if (n == zmath::nzinf) return "0"; // 0
        if (n == 0) return std::to_string(coef_[0]); // 常数

        // 是不是负数
        auto is_negative = [](const double& d) -> bool {
            if (d < 0) return true;
            return false;
        };

        std::string str = "";
        std::string leader = std::to_string(coef_[0]) + " x^" + std::to_string(n); // 首项
        str += leader;
        for (int i = 1; i <= n - 1; i++) {
            double current_coef = coef_[i];

            if (current_coef == 0)
                continue;

            if (is_negative(current_coef)) {
                str += " - "; // 先加上负号
                current_coef = -current_coef; // 取反 
                // 这样做是为了好看一点
            } else {
                str += " + ";
            }
            
            if (i != n - 1) {
                str += std::to_string(current_coef) + " x^" + std::to_string(n - i);
            } else {
                str += std::to_string(current_coef) + " x";
            }
        }

        double last_coef = coef_.back();
        if (last_coef != 0) {
            if (is_negative(last_coef)) {
                str += " - ";
                str += std::to_string(-last_coef);
            } else {
                str += " + ";
                str += std::to_string(last_coef);
            }
        }

        return str;
    }

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& poly) {
        os << poly.to_string();
        return os;
    }

    void print() const {
        std::cout << *this << std::endl;
    }

    // TODO: 最小二乘法需要解线性方程组, 等linalg写完再回来补
    /**
     * @brief 最小二乘法拟合多项式
     * @param x 横坐标
     * @param y 纵坐标
     * @param order 拟合的多项式的度
     * @return Polynomial 
     */
    static Polynomial fit(const std::vector<double>& x, const std::vector<double>& y, int order);

private:
    std::vector<double> coef_;
};

}