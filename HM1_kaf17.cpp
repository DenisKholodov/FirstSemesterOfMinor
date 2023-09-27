#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>

#define _USE_MATH_DEFINES
#include <cmath>

void SineSignal(int, double, double, double = 10., int = 1024);
void TriangularSignal (int, double, double, double, double = 10., int = 1024);
void RectaingularSignal (int, double, double, double, double, double = 10., int = 1024);
void TwoRandomSignal (int, double = 10., int = 1024);

int main()
{
    SineSignal(60, 0., 2*M_PI);
    TriangularSignal(100, 3., 0., 4);
    RectaingularSignal(100, 5.5, 2, 0., 3);
    TwoRandomSignal(41);
    return 0;
}

void SineSignal (int N_samples, double phi0, double T, double F, int quantization_levels_num)
{
    double quantization_min = -1., quantization_max = 1.;
    const char csv_file_name[64] = "data_sine.csv"; // имя файла для вывода

    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n"; // записать заголовки колонок

    for (int i = 0; i < N_samples; ++i)
    {
        double signal_val = trunc((sin((i/F)*2*M_PI/T + phi0)*(quantization_levels_num - 1))/(quantization_max - quantization_min)) * \
            ((quantization_max - quantization_min)/(quantization_levels_num - 1));
        csv_file << (i / F) << "," << signal_val << "\n";
    }
    csv_file.close();
}

void TriangularSignal (int N_samples, double tau, double phi0, double T, double F, int quantization_levels_num)
{
    double quantization_min = 0., quantization_max = tau/2;
    const char csv_file_name[64] = "data_triangal.csv"; // имя файла для вывода
    double digital_signal[N_samples] = {};

    int k = 0;
    for (int i = 0; i < N_samples; i++)
    {
        double t = i/F - k*T;

        if (tau/2 + phi0 - t > 0.0000001){
            digital_signal[i] = t;
        }
        else if (tau + phi0 - t > 0.0000001){
            digital_signal[i] = 2*quantization_max - t;
        }
        else
        {
            digital_signal[i] = 0.;
        }

        if (i/F - T*(k+1) > 0.0000001) k++;
    }

    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n"; // записать заголовки колонок

    for (int i = 0; i < N_samples; ++i)
    {
        double signal_val = trunc((digital_signal[i]*(quantization_levels_num - 1))/(quantization_max - quantization_min)) * \
            ((quantization_max - quantization_min)/(quantization_levels_num - 1));
        csv_file << (i / F) << "," << signal_val << "\n";
    }
    csv_file.close();
}

void RectaingularSignal (int N_samples, double MaxFunc, double tau, double phi0, double T, double F, int quantization_levels_num)
{
    double quantization_min = 0., quantization_max = MaxFunc;
    const char csv_file_name[64] = "data_rectangal.csv"; // имя файла для вывода
    double digital_signal[N_samples] = {};

    int k = 0;
    for (int i = 0; i < N_samples; i++)
    {
        double t = i/F - k*T;

        if (tau + phi0 - t > 0.0000001){
            digital_signal[i] = MaxFunc;
        }
        else
        {
            digital_signal[i] = 0.;
        }

        if (i/F - T*(k+1) > 0.0000001) k++;
    }
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n"; // записать заголовки колонок

    for (int i = 0; i < N_samples; ++i)
    {
        double signal_val = trunc((digital_signal[i]*(quantization_levels_num - 1))/(quantization_max - quantization_min)) * \
            ((quantization_max - quantization_min)/(quantization_levels_num - 1));
        csv_file << (i / F) << "," << signal_val << "\n";
    }
    csv_file.close();
}

void TwoRandomSignal (int N_samples, double F, int quantization_levels_num)
{
    double quantization_min = 0., quantization_max = 1.;
    const char csv_file_name[64] = "data_randomTwo.csv"; // имя файла для вывода
    double digital_signal[N_samples] = {};

    srand (static_cast <unsigned> (time(0)));

    double X = 1.;
    
    double r1 = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/X));
    double r2 = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/X)) + r1;
    double r3 = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/X)) + r2;
    double r4 = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/X)) + r3;
    for (int i = 0; i < N_samples; i++)
    {
        double t = i/F;
        if (r1 - t > 0.0000001){
            digital_signal[i] = 0;
        }
        else if (r2 - t > 0.0000001)
        {
            digital_signal[i] = quantization_max;
        }
        else if (r3 - t > 0.0000001)
        {
            digital_signal[i] = 0;
        }
        else if (r4 - t > 0.0000001)
        {
            digital_signal[i] = quantization_max;
        }
        else
        {
            digital_signal[i] = 0;
        }
    }
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal\n"; // записать заголовки колонок

    for (int i = 0; i < N_samples; ++i)
    {
        double signal_val = trunc((digital_signal[i]*(quantization_levels_num - 1))/(quantization_max - quantization_min)) * \
            ((quantization_max - quantization_min)/(quantization_levels_num - 1));
        csv_file << (i / F) << "," << signal_val << "\n";
    }
    csv_file.close();
    std::cout << r1 << '\n' << r2 << '\n' << r3 << '\n' << r4 << '\n';
    }