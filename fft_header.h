#pragma once

#include <complex>
#include <vector>

class FastFourierTransform {
public:
    using ComplexType = std::complex<double>;

    // Основной метод для вычисления БПФ
    static std::vector<ComplexType> calculate(std::vector<ComplexType>& input, bool inverseFlag);

private:
    // Вспомогательные методы
    static int reverse_bits(int value, int bit_count); // Битовая перестановка для индексирования

    // Универсальный метод FFT (например, Bluestein или Cooley-Tukey)
    static std::vector<ComplexType> fft_mixed(const std::vector<ComplexType>& input, bool inverseFlag);

    // Специализированные реализации FFT для определённых простых чисел
    static std::vector<ComplexType> fft_2(const std::vector<ComplexType>& input, bool inverseFlag);
    static std::vector<ComplexType> fft_3(const std::vector<ComplexType>& input, bool inverseFlag);
    static std::vector<ComplexType> fft_5(const std::vector<ComplexType>& input, bool inverseFlag);

    // Прямое вычисление DFT для малых последовательностей
    static std::vector<ComplexType> dft(const std::vector<ComplexType>& input, bool inverseFlag);
};