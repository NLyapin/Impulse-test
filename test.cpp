#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>

using namespace std;

// --- ВАРИАНТ 2 ---

class FFT {
public:
    FFT(int N) : N(N) {
        // Проверка, что размер кратен 2, 3 или 5
        if (N <= 0 || (N % 2 != 0 && N % 3 != 0 && N % 5 != 0)) {
            throw invalid_argument("Размер должен быть кратен 2, 3 или 5.");
        }
    }

    // Прямое преобразование Фурье
    void forward(vector<complex<double>>& data) {
        fft(data, false);
    }

    // Обратное преобразование Фурье
    void inverse(vector<complex<double>>& data) {
        fft(data, true);
    }

private:
    int N;

    // Основная функция для вычисления FFT
    void fft(vector<complex<double>>& data, bool inverse) {
        int n = data.size();
        if (n <= 1) return;

        vector<complex<double>> even(n / 2), odd(n / 2);

        // Разделение на чётные и нечётные индексы
        for (int i = 0; i < n / 2; i++) {
            even[i] = data[2 * i];
            odd[i] = data[2 * i + 1];
        }

        // Рекурсивный вызов FFT
        fft(even, inverse);
        fft(odd, inverse);

        double angle = 2 * M_PI / n * (inverse ? -1 : 1);
        complex<double> w(1), wn(cos(angle), sin(angle));
        
        // Вычисление финальных значений
        for (int i = 0; i < n / 2; i++) {
            data[i] = even[i] + w * odd[i];
            data[i + n / 2] = even[i] - w * odd[i];
            if (inverse) {
                data[i] /= 2;
                data[i + n / 2] /= 2;
            }
            w *= wn;
        }
    }
};

// Функция для генерации случайных данных
vector<complex<double>> generate_random_data(int N) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(-1.0, 1.0);
    vector<complex<double>> data(N);

    for (int i = 0; i < N; i++) {
        data[i] = complex<double>(dis(gen), dis(gen));
    }
    return data;
}

// Функция для вычисления ошибки между оригинальными и преобразованными данными
double compute_error(const vector<complex<double>>& original, const vector<complex<double>>& transformed) {
    double error = 0.0;
    for (size_t i = 0; i < original.size(); i++) {
        error += abs(original[i] - transformed[i]);
    }
    return error / original.size();
}

int main() {
    try {
        int N = 8; // Размер преобразования
        FFT fft(N);
        
        // Генерация случайных данных
        vector<complex<double>> data = generate_random_data(N);
        vector<complex<double>> original_data = data;

        // Прямое и обратное преобразования
        fft.forward(data);
        fft.inverse(data);

        // Вычисление ошибки
        double error = compute_error(original_data, data);
        
        cout << "Ошибка между входными и выходными данными: " << error << endl;
    }
    catch (const exception& e) {
        cout << "Ошибка: " << e.what() << endl;
    }

    return 0;
}