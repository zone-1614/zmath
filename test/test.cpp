// 不知道那些自动化测试是啥, 自己写点测试的代�?
#include <iostream>
#include <zmath.h>

using namespace zmath;

complex normalize(complex c) {
    if (abs(c.real()) < epsilon) {
        c.real(0);
    }
    if (abs(c.imag() < epsilon)) {
        c.imag(0);
    }
    return c;
}

void test_fft() {
    std::cout << ninf << std::endl;
    complex c[] = { 3, 2, 1, 0 };
    complex_array coef(c, 4);
    fft(coef);
    for (int i = 0; i < 4; i++) {
        std::cout << normalize(coef[i]) << std::endl;
    }
    ifft(coef);
    for (int i = 0; i < 4; i++) {
        std::cout << normalize(coef[i]) << std::endl;
    }
}

void test_polynomial() {
    // 什么都没有的多项式, 默认�? 0
    std::vector<double> coef0 { };
    Polynomial poly0(coef0);
    poly0.print();

    // 默认构�?
    std::vector<double> coef1;
    Polynomial poly1(coef1);
    poly1.print();

    std::vector<double> coef2{ 0, 0, 0, 1, 2, 3, 0, 0, 1 };
    Polynomial poly2(coef2);
    poly2.print();

    // 测试 print
    std::cout << std::endl << "测试 print " << std::endl;
    std::vector<double> coef { -3.5, 0, -5.3, 2.2, -6.5, -9.999 }; 
    // -3.500000 x^5 - 5.300000 x^3 + 2.200000 x^2 - 6.500000 x - 9.999000
    Polynomial poly(coef);
    poly.print();

    // 求导
    std::cout << std::endl << "测试求导 derivative " << std::endl;
    auto deriv = poly.derivative();
    deriv.print();

    // 首一
    std::cout << std::endl << "测试首一多项�? monic " << std::endl;
    auto monic = deriv.monic();
    monic.print();
    // 测试设置系数
    std::cout << std::endl << "测试设置系数" << std::endl;
    monic.set_coef(3, 2.6);
    monic.print();
    monic.set_coef(0, 3.33);
    monic.set_coef(3, 0);
    monic.set_coef(4, 0);
    monic.print();

    // 测试 operator()
    std::cout << std::endl << "测试 operator()" << std::endl;
    std::vector<double> coef3 { 1, -1, 2 }; 
    Polynomial poly3(coef3);
    poly3.print();
    std::cout << poly3(0) << " " << poly3(1) << " " << poly3(2) << std::endl;
    std::cout << poly3(-1) << " " << poly3(3) << " " << poly3(0.5) << std::endl;

    // 测试 operator+ operator-
    std::cout << std::endl << "测试 operator+ operator-" << std::endl;
    std::vector<double> coef_add1{ 1, 2, 3, 0, -6, 1 };
    Polynomial poly_add1(coef_add1);
    poly_add1.print();
    std::vector<double> coef_add2{ -1, -2, 3, 4, 0, 1 };
    Polynomial poly_add2(coef_add2);
    poly_add2.print();

    auto poly_add = poly_add1 + poly_add2;
    poly_add.print();

    auto poly_subtract = poly_add1 - poly_add2;
    poly_subtract.print();

    // 测试 operator* (数乘)
    auto mul1 = poly_add * 2;
    mul1.print();
    auto mul2 = -2.1 * poly_add;
    mul2.print();

    // 测试多项式乘法
    std::cout << std::endl << "测试多项式乘法" << std::endl;
    std::vector<double> coef_mul1{ 2, -4, 0.5, -1 };
    Polynomial poly_mul1(coef_mul1);
    poly_mul1.print();
    std::vector<double> coef_mul2{ -1, 0, 3 };
    Polynomial poly_mul2(coef_mul2);
    poly_mul2.print();
    auto poly_mul = poly_mul1 * poly_mul2;
    poly_mul.print();

    // 测试多项式幂
    std::cout << std::endl << "测试多项式幂" << std::endl;
    auto poly_pow = poly_mul2^0;
    auto poly_pow1 = poly_mul2^1;
    auto poly_pow2 = poly_mul2^2;
    poly_pow.print();
    poly_pow1.print();
    poly_pow2.print();

    // 测试 += -= *= ^=
    std::cout << std::endl << "测试 += -= *= ^=" << std::endl;
    std::vector<double> coef_z1{ 2, -4, 0.5, -1 };
    Polynomial poly_z1(coef_z1);
    poly_z1.print();
    std::vector<double> coef_z2{ -1, 0, 3 };
    Polynomial poly_z2(coef_z2);
    poly_z2.print();
    poly_z1 += poly_z2;
    poly_z1.print();
    poly_z1 *= -2.5;
    poly_z1.print();
    poly_z1 *= poly_z2;
    poly_z1.print();
    poly_z2^=2;
    poly_z2.print();
}

void test_linalg() {
    // 构造并打印
    std::vector<double> vv{ 1, 2, 3 };
    Vector v(vv);
    v.print();
    
    // 转置
    auto vt = v.transpose();
    vt.print();

    auto va = vt.transpose() + v;
    va.print();

    auto vb = va - 2.5 * v;
    vb.print();
}

void test_vec2() {
    Vec2 v1(2.0, 2.0);
    Vec2 v2(1.0, -1.0);
    v1.print();
    fmt::print("v1 v2 angle {:.3f}\nv1 norm {:.3f}\nv1 v2 distance {:.3f}\n", v1.angle(v2), v1.norm(), v1.distance(v2));
    auto v = v1.normalize();
    v.print();
}

int main() {
    // test_fft();
    // test_polynomial();
    // test_linalg();
    test_vec2();
}