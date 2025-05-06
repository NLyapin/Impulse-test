#include "fft_header.h"
#include <cmath>
#include <stdexcept>
#include <functional>

constexpr double kPi = 3.14159265358979323846;

// Основной метод для расчёта FFT
std::vector<FastFourierTransform::ComplexType> FastFourierTransform::calculate(std::vector<ComplexType>& input, bool inverseFlag) {
    auto size = input.size();
    auto transformed = fft_mixed(input, inverseFlag);

    if (inverseFlag) {
        // Для обратного преобразования делим результаты на размер входных данных
        for (auto& elem : transformed) {
            elem /= static_cast<double>(size);
        }
    }

    return transformed;
}

// Функция для инвертирования битов
int FastFourierTransform::reverse_bits(int value, int bit_count) {
    int result = 0;
    for (int i = 0; i < bit_count; ++i) {
        result = (result << 1) | ((value >> i) & 1);
    }
    return result;
}

// Метод для выбора оптимального алгоритма FFT в зависимости от размера входных данных
std::vector<FastFourierTransform::ComplexType> FastFourierTransform::fft_mixed(const std::vector<ComplexType>& input, bool inverseFlag) {
    int n = input.size();
    if (n % 2 == 0) {
        return fft_2(input, inverseFlag); // FFT для размера, кратного 2
    } else if (n % 3 == 0) {
        return fft_3(input, inverseFlag); // FFT для размера, кратного 3
    } else if (n % 5 == 0) {
        return fft_5(input, inverseFlag); // FFT для размера, кратного 5
    } else {
        return dft(input, inverseFlag); // Обыкновенное дискретное преобразование Фурье
    }
}

// Алгоритм FFT для размера, кратного 2
std::vector<FastFourierTransform::ComplexType> FastFourierTransform::fft_2(const std::vector<ComplexType>& input, bool inverseFlag) {
    int n = input.size();
    std::vector<ComplexType> even(n / 2), odd(n / 2);

    // Разделение на чётные и нечётные элементы
    for (int i = 0; i < n / 2; ++i) {
        even[i] = input[2 * i];
        odd[i] = input[2 * i + 1];
    }

    auto fft_even = fft_mixed(even, inverseFlag);
    auto fft_odd = fft_mixed(odd, inverseFlag);

    std::vector<ComplexType> output(n);
    double sign = inverseFlag ? 1.0 : -1.0;

    // Комбинирование результатов чётной и нечётной части
    for (int k = 0; k < n / 2; ++k) {
        double angle = 2.0 * kPi * sign * k / n;
        ComplexType twiddle = std::polar(1.0, angle);
        output[k] = fft_even[k] + twiddle * fft_odd[k];
        output[k + n / 2] = fft_even[k] - twiddle * fft_odd[k];
    }

    return output;
}

// Алгоритм FFT для размера, кратного 3
std::vector<FastFourierTransform::ComplexType> FastFourierTransform::fft_3(const std::vector<ComplexType>& input, bool inverseFlag) {
    int n = input.size();
    int step = n / 3;

    std::vector<ComplexType> a(step), b(step), c(step);
    for (int i = 0; i < step; ++i) {
        a[i] = input[3 * i];
        b[i] = input[3 * i + 1];
        c[i] = input[3 * i + 2];
    }

    auto A = fft_mixed(a, inverseFlag);
    auto B = fft_mixed(b, inverseFlag);
    auto C = fft_mixed(c, inverseFlag);

    std::vector<ComplexType> result(n);
    int direction = inverseFlag ? 1 : -1;

    // Комбинирование результатов для трёх частей
    for (int k = 0; k < step; ++k) {
        for (int i = 0; i < 3; ++i) {
            int shift = k + i * step;
            ComplexType w1 = std::polar(1.0, 2.0 * kPi * direction * shift / n);
            ComplexType w2 = std::polar(1.0, 4.0 * kPi * direction * shift / n);
            result[shift] = A[k] + w1 * B[k] + w2 * C[k];
        }
    }

    return result;
}

// Алгоритм FFT для размера, кратного 5
std::vector<FastFourierTransform::ComplexType> FastFourierTransform::fft_5(const std::vector<ComplexType>& input, bool inverseFlag) {
    int n = input.size();
    int step = n / 5;

    std::vector<ComplexType> a(step), b(step), c(step), d(step), e(step);
    for (int i = 0; i < step; ++i) {
        a[i] = input[5 * i];
        b[i] = input[5 * i + 1];
        c[i] = input[5 * i + 2];
        d[i] = input[5 * i + 3];
        e[i] = input[5 * i + 4];
    }

    auto A = fft_mixed(a, inverseFlag);
    auto B = fft_mixed(b, inverseFlag);
    auto C = fft_mixed(c, inverseFlag);
    auto D = fft_mixed(d, inverseFlag);
    auto E = fft_mixed(e, inverseFlag);

    std::vector<ComplexType> result(n);
    int direction = inverseFlag ? 1 : -1;

    // Комбинирование результатов для пяти частей
    for (int k = 0; k < step; ++k) {
        for (int i = 0; i < 5; ++i) {
            int index = k + i * step;
            ComplexType sum = A[k];
            for (int j = 1; j < 5; ++j) {
                double angle = 2.0 * kPi * direction * j * index / n;
                sum += std::polar(1.0, angle) * (j == 1 ? B[k] : j == 2 ? C[k] : j == 3 ? D[k] : E[k]);
            }
            result[index] = sum;
        }
    }

    return result;
}

// Обыкновенное дискретное преобразование Фурье (DFT)
std::vector<FastFourierTransform::ComplexType> FastFourierTransform::dft(const std::vector<ComplexType>& input, bool inverseFlag) {
    int n = input.size();
    std::vector<ComplexType> output(n);
    double direction = inverseFlag ? 1.0 : -1.0;

    // Стандартное вычисление DFT
    for (int i = 0; i < n; ++i) {
        ComplexType sum = 0.0;
        for (int j = 0; j < n; ++j) {
            double angle = 2.0 * kPi * direction * i * j / n;
            sum += input[j] * std::polar(1.0, angle);
        }
        output[i] = sum;
    }

    return output;
}
