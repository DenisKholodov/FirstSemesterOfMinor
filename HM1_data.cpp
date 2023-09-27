#include <fstream>
#include <iostream>
#include <math.h>
#include <cstdint>

int main()
{
    double F = 10.; // частота дискретизации
    int quantization_levels_num = 1024; // количество уровней равномерного квантования
    double quantization_min = -1., quantization_max = 1.;
    int N_samples = 1024;
    int digital_signal[1024] = {}; // коды квантования цифрового сигнала
    const char csv_file_name[64] = "data.csv"; // имя файла для вывода

    for (int i = 0; i < N_samples; ++i)
    {
        digital_signal[i] = i;
    }

    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n"; // записать заголовки колонок
    for (int i = 0; i < N_samples; ++i)
    {
        // перевод кода равномерного квантования в физическую величину
        double signal_val = ((quantization_max - quantization_min) * \
            double(digital_signal[i]) / (quantization_levels_num - 1)) + quantization_min;
        // вывод через запятую значений колонок
        csv_file << (i / F) << "," << signal_val << "\n";
    }
    csv_file.close();

    return 0;
}