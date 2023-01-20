#include <iostream>
#include <zmath/fft.h>

zmath::complex normalize(zmath::complex c) {
    if (abs(c.real()) < zmath::epsilon) {
        c.real(0);
    }
    if (abs(c.imag() < zmath::epsilon)) {
        c.imag(0);
    }
    return c;
}

void test_fft() {
    zmath::complex c[] = { 3, 2, 1, 0 };
    zmath::complex_array coef(c, 4);
    zmath::fft(coef);
    for (int i = 0; i < 4; i++) {
        std::cout << normalize(coef[i]) << std::endl;
    }
    zmath::ifft(coef);
    for (int i = 0; i < 4; i++) {
        std::cout << normalize(coef[i]) << std::endl;
    }
}

int main() {
    test_fft();
}