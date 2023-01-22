#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <stdexcept>
#include <memory>
#include <variant>

namespace zmath {

enum class VecType {
    Row = 0,
    Col = 1
};

// 乘法的返回类型
using MulResult = std::variant<double, std::shared_ptr<Vector>, std::shared_ptr<Matrix>>; 

class Matrix;

class Vector {
public:
    friend class Matrix;

    Vector(size_t n = 1, VecType type = VecType::Col)
        : data_(n, 0), type_(type) { }

    Vector(const std::vector<double>& data, VecType type = VecType::Col) 
        : data_(data), type_(type) { }

    VecType type() const {
        return type_;
    }

    void set_type(VecType type) {
        type_ = type;
    }

    size_t size() const {
        return data_.size();
    }

    // 转置
    Vector transpose() const {
        Vector ret = *this;
        if (ret.type_ == VecType::Col) {
            ret.type_ = VecType::Row;
        } else {
            ret.type_ = VecType::Col;
        }
        return ret;
    }

    double& operator[](size_t i) {
        return data_[i];
    }

    const double& operator[](size_t i) const {
        return data_[i];
    }

    Vector operator+(const Vector& rhs) const {
        check_type(*this, rhs);
        check_size(*this, rhs);

        Vector ret = *this;
        for (int i = 0; i < ret.size(); i++) {
            ret[i] += rhs[i];
        }
        return ret;
    }

    Vector operator-(const Vector& rhs) const {
        check_type(*this, rhs);
        check_size(*this, rhs);

        Vector ret = *this;
        for (int i = 0; i < ret.size(); i++) {
            ret[i] -= rhs[i];
        }
        return ret;
    }

    Vector operator*(double d) const {
        Vector ret = *this;
        for (int i = 0; i < ret.size(); i++) {
            ret[i] *= d;
        }
        return ret;
    }

    friend Vector operator*(double d, const Vector& rhs) {
        Vector ret = rhs;
        for (int i = 0; i < ret.size(); i++) {
            ret[i] *= d;
        }
        return ret;
    }

    // 乘法需要按照不同类型产生不同结果, 可以是内积, 数, 矩阵
    // Vector operator*(const Vector& rhs) const {
    //     check_size(*this, rhs);
        
    //     Vector ret = *this;
    //     for (int i = 0; i < ret.size(); i++) {
    //         ret[i] *= rhs[i];
    //     }
    //     return ret;
    // }

    std::string to_string() const {
        std::string ret = "";
        if (type_ == VecType::Col)
            ret += "Col ( ";
        else 
            ret += "Row ( ";
        
        for (int i = 0; i < data_.size(); i++) {
            ret += std::to_string(data_[i]);
            if (i == data_.size() - 1)
                ret += " )";
            else 
                ret += ", ";
        }
        return ret;
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector& v) {
        os << v.to_string();
        return os;
    }

    void print() const {
        std::cout << *this << std::endl;
    }

private: 
    std::vector<double> data_;
    VecType type_;

    static void check_size(const Vector& lhs, const Vector& rhs) {
        if (lhs.size() != rhs.size())
            throw std::logic_error("vector size not equal");
    }

    static void check_type(const Vector& lhs, const Vector& rhs) {
        if (lhs.type() != rhs.type())
            throw std::logic_error("vector type not equal");
    }
};

}