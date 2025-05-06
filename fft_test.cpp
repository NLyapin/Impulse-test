#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <iomanip>

#include "fft_header.h"

int main() {
    using Cmplx = std::complex<double>;
    const int length = 30;  // длина сигнала, кратная 2, 3 и 5

    std::vector<Cmplx> signal(length);
    std::mt19937 gen(42);  // инициализация генератора с фиксированным зерном
    std::uniform_int_distribution<int> binary(0, 1);

    // Генерация QPSK: 2 бита -> 1 символ
    auto encode_qpsk = [](int high, int low) -> Cmplx {
        const double scale = std::sqrt(0.5);  // нормировка
        double re = high ? -1.0 : 1.0;
        double im = low ? -1.0 : 1.0;
        return Cmplx(re, im) * scale;
    };

    for (int i = 0; i < length; ++i) {
        int bit_hi = binary(gen);
        int bit_lo = binary(gen);
        signal[i] = encode_qpsk(bit_hi, bit_lo);
    }

    std::cout << "\nРазмер сигнала: " << length << "\n";

    std::cout << "\nСимволы QPSK:\n";
    for (const auto& sym : signal) {
        std::cout << std::fixed << std::setprecision(3)
                  << "(" << sym.real() << ", " << sym.imag() << "i)\n";
    }

    // Прямое БПФ
    std::vector<Cmplx> spectrum = FastFourierTransform::calculate(signal, false);

    std::cout << "\nСпектр (FFT):\n";
    for (const auto& freq : spectrum) {
        std::cout << std::fixed << std::setprecision(3)
                  << "(" << freq.real() << ", " << freq.imag() << "i)\n";
    }

    // Обратное БПФ
    std::vector<Cmplx> restored = FastFourierTransform::calculate(spectrum, true);

    std::cout << "\nВосстановленный сигнал (IFFT):\n";
    for (const auto& sample : restored) {
        std::cout << std::fixed << std::setprecision(3)
                  << "(" << sample.real() << ", " << sample.imag() << "i)\n";
    }

    // Расчет среднеквадратичной ошибки
    double mse = 0.0;
    for (int i = 0; i < length; ++i) {
        double real_diff = signal[i].real() - restored[i].real();
        double imag_diff = signal[i].imag() - restored[i].imag();
        mse += real_diff * real_diff + imag_diff * imag_diff;
    }
    mse /= length;

    std::cout << "\nСреднеквадратичная ошибка (MSE): " << std::scientific << mse << "\n";

    return 0;
}
