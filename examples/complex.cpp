// Complex numbers and polynomial roots with MathRix.

#include <iostream>
#include <mathrix/mathrix.hpp>

using namespace mathrix;
using namespace mathrix::literals;

int main() {
    const Complex z1 = 3.0 + 4.0_i;
    const Complex z2{1, -2};
    std::cout << "z1 = " << z1 << ", z2 = " << z2 << "\n";
    std::cout << "z1 + z2 = " << z1 + z2 << "\n";
    std::cout << "z1 - z2 = " << z1 - z2 << "\n";
    std::cout << "z1 * z2 = " << z1 * z2 << "\n";
    std::cout << "z1 / z2 = " << z1 / z2 << "\n";
    std::cout << "|z1| = " << z1.abs() << ", arg z1 = " << z1.arg() << " rad, conj z1 = " << z1.conj() << "\n";
    std::cout << "e^(i pi) = " << exp(Complex(0, kPi)) << "  (Euler)\n";
    std::cout << "sqrt(-4) = " << sqrt(Complex(-4)) << ", i^i = " << pow(1.0_i, 1.0_i) << "\n";

    std::cout << "\nfourth roots of 16:";
    for (const Complex& w : nthRoots(Complex(16), 4)) std::cout << "  " << chop(w);

    // x^3 - 2x^2 + 4x - 8 = (x - 2)(x^2 + 4)  ->  roots 2, 2i, -2i
    const Polynomial p{-8, 4, -2, 1};
    std::cout << "\n\nroots of " << p << ":";
    for (const Complex& r : p.roots()) std::cout << "  " << r;
    std::cout << "\n";
}
