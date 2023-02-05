#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <variant>

#include <fmt/format.h>

namespace zmath {

enum class VecType {
    Row = 0,
    Col = 1
};


class Vector;
class Matrix;

// 乘法的返回类型
using MulResult = std::variant<double, std::shared_ptr<Matrix>>; 

class Vector {
public:
    friend class Matrix;
    friend MulResult;

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

    double& operator()(size_t i) {
        return data_[i];
    }

    const double& operator()(size_t i) const {
        return data_[i];
    }

    Vector& operator+=(const Vector& rhs) {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] += rhs.data_[i];
        }
        return *this;
    }

    Vector& operator-=(const Vector& rhs) {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] -= rhs.data_[i];
        }
        return *this;
    }

    Vector& operator*=(double c) {
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] *= c;
        }
        return *this;
    }

    friend Vector operator*(Vector lhs, double c) {
        lhs *= c;
        return lhs;
    }

    friend Vector operator*(double c, Vector rhs) {
        rhs *= c;
        return rhs;
    }

    friend Vector operator+(Vector lhs, const Vector& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend Vector operator-(Vector lhs, const Vector& rhs) {
        lhs -= rhs;
        return lhs;
    }

    std::string to_string() const {
        return fmt::format("{} ( {:.3f} )", type_ == VecType::Col ? "Col" : "Row", fmt::join(data_, ", "));
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
};


class Matrix {
public:
    friend MulResult;
    friend MulResult operator*(const Vector& lhs, const Vector& rhs);

    using PLU_Type = std::pair<std::vector<size_t>, Matrix>;

    Matrix(size_t row = 1, size_t col = 1);
    Matrix(const std::vector<std::vector<double>>& data);

    bool is_squared() const {
        return get_row_size() == get_col_size();
    }
    bool is_symmetric() const {
        if (!is_squared()) return false;

        auto n = get_row_size();
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                if (!eq(data_[i][j], data_[j][i]))
                    return false;
            }
        }
        return true;
    }

    size_t get_row_size() const {
        return data_.size();
    }
    size_t get_col_size() const {
        return data_[0].size();
    }

    Vector get_n_row_vector(size_t n) const {
        return Vector(data_[n], VecType::Row);
    }
    Vector get_n_col_vector(size_t n) const {
        auto s = get_row_size();
        Vector v(s);
        for (int i = 0; i < s; i++) {
            v(i) = data_[i][n];
        }
        return v;
    }

    Matrix LU_decomp() const {
        // TODO if (!is_squared()) error
        Matrix LU(*this);

        auto n = get_row_size();
        for (int j = 0; j < n - 1; j++) {
            // L
            double cj = 1.0 / LU(j, j);
            for (int i = j + 1; i < n; i++) {
                LU(i, j) *= cj;
            }
            // U
            for (int i = j + 1; i < n; i++) {
                for (int k = j + 1; k < n; k++) {
                    LU(i, k) -= LU(i, j) * LU(j, k);
                }
            }
        }
        return LU;
    }

    Vector LU_solve(const Vector& b) const {
        // TODO if (!is_squared()) error
        Vector x(b.size());
        Matrix LU = LU_decomp();
        auto y = LU_solve_L(LU, b);
        x = LU_solve_U(LU, y);
        return x;
    }

    double det() const {
        // TODO if (!is_squared()) error
        Matrix LU = LU_decomp();
        double det = 1.0;
        for (int i = 0; i < get_row_size(); i++) {
            det *= LU(i, i);
        }
        return det;
    }

    Matrix inv() const {
        // TODO if (!is_squared()) error
        auto n = get_row_size();
        Matrix mat(n, n);
        Matrix LU = LU_decomp();
        Vector ej(n), xj(n);
        ej(0) = 1.0;
        for (int j = 0; j < n; j++) {
            if (j >= 1) std::swap(ej(j-1), ej(j));
            auto yj = LU_solve_L(LU, ej);
            xj = LU_solve_U(LU, yj);
            for (int i = 0; i < n; i++) {
                mat(i, j) = xj(i);
            }
        }
        return mat;
    }

    Matrix T() const {
        const int r = get_row_size(), c = get_col_size();
        Matrix mat(c, r);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                mat(j, i) = data_[i][j];
            }
        }
        return mat;
    }

    const double& operator()(size_t i, size_t j) const;
    double& operator()(size_t i, size_t j);
    Matrix& operator+=(const Matrix& rhs);
    Matrix& operator-=(const Matrix& rhs); 
    Matrix& operator*=(double c);

    friend Matrix operator*(Matrix lhs, double c);
    friend Matrix operator*(double c, Matrix rhs);
    friend Matrix operator+(Matrix lhs, const Matrix& rhs);
    friend Matrix operator-(Matrix lhs, const Matrix& rhs);
    friend Matrix operator*(Matrix lhs, const Matrix& rhs);

    std::string to_string() const {
        std::string str = "Mat [ ";
        const size_t m = get_row_size();
        for (int i = 0; i < m; i++) {
            str += fmt::format("{:.3f}", fmt::join(data_[i], ", "));
            if (i != m - 1) {
                str += ";\n ";
            }
        }
        str += " ]";
        return str;
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix& m) {
        os << m.to_string();
        return os;
    }

    void print() const {
        std::cout << *this << std::endl;
    }

private: 
    std::vector<std::vector<double>> data_;

    Vector LU_solve_L(const Matrix& LU, const Vector& b) const {
        auto n = b.size();
        Vector x(n);
        for (int i = 0; i < n; i++) {
            x(i) = b(i);
            for (int j = 0; j < i; j++) {
                x(i) -= LU(i, j) * x(j);
            }
        }
        return x;
    }

    Vector LU_solve_U(const Matrix& LU, const Vector& b) const {
        auto n = b.size();
        Vector x(n);
        for (int i = n - 1; i >= 0; i--) {
            x(i) = b(i);
            for (int j = i + 1; j < n; j++) {
                x(i) -= LU(i, j) * x(j);
            }
            x(i) /= LU(i, i);
        }
        return x;
    }
};

MulResult operator*(const Vector& lhs, const Vector& rhs) {
    MulResult result;
    if (lhs.type() == VecType::Row && rhs.type() == VecType::Col) {
        double r = 0.0;
        for (int i = 0; i < lhs.size(); i++) {
            r += lhs(i) * rhs(i);
        }
        result = r;
        return result;
    } else if (lhs.type() == VecType::Col && rhs.type() == VecType::Row) {
        auto mat = std::make_shared<Matrix>(lhs.size(), rhs.size());
        for (int i = 0; i < lhs.size(); i++) {
            for (int j = 0; j < rhs.size(); j++) {
                mat->data_[i][j] = lhs(i) * rhs(j);
            }
        }
        result = mat;
        return result;
    } else {
        result = zmath::inf; // error
        return result;
    }
}

}