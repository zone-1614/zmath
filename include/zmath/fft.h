#pragma once
#include <complex>
#include <valarray>
#include <zmath/constant.h>

// fft 用于计算多项式乘积

namespace zmath {

using complex = std::complex<double>;
using complex_array = std::valarray<complex>;
using zmath::pi;

/**
 * @brief 
 * 
 * @param coef 多项式的系数. 例如, 如果多项式是 1 - 2x + 3x^2, 则 coef 是 { 1, -2, 3 } 
 */
void fft(complex_array& coef) {
    const int N = coef.size(); // degree of 

    if (N <= 1) return;

    complex_array even = coef[std::slice(0, N / 2, 2)];
    complex_array odd  = coef[std::slice(1, N / 2, 2)];

    fft(even);
    fft(odd);

    for (size_t k = 0; k < N / 2; k++) {
        complex t = std::polar(1.0, -2 * pi * k / N) * odd[k];
        coef[k] = even[k] + t;
        coef[k + N / 2] = even[k] - t;
    }
}

void ifft(complex_array& coef) {
    coef = coef.apply(std::conj);
    fft(coef);
    coef = coef.apply(std::conj);
    coef /= coef.size();
}

}