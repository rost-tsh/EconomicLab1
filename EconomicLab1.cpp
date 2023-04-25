// EconomicLab1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <random>

using namespace std;

const double PI = 3.141592653589793;
const double a = -3, b = 3;
const int N = 1000;
const int M = 10000;

double f(double x)
{
    return sin(x);
}

double p(double x)
{
    return 1/sqrt(2 * PI) * exp(-pow(x,2)/2);
}

double CentralIntegral(double a, double b, int n)
{
    double Integral = 0.0;
    double h = (b - a) / n;
    for (int i = 1; i <= n; i++)
        Integral = Integral + h * f(a + h * (i - 0.5));
    
}

double TrapezeIntegral(double a, double b, int n)
{
    double h = (b - a) / n;
    double Integral = h * (f(a) + f(b)) / 2.0;
    for (int i = 1; i <= n - 1; i++)
        Integral = Integral + h * f(a + h * i);
    
}

double SympsonIntegral(double a, double b, int n)
{
    double h = (b - a) / n;
    double Integral = h * (f(a) + f(b)) / 6.0;
    for (int i = 1; i <= n; i++)
        Integral = Integral + 4.0 / 6.0 * h * f(a + h * (i - 0.5));
    for (int i = 1; i <= n - 1; i++)
        Integral = Integral + 2.0 / 6.0 * h * f(a + h * i);
    
}

double p(double x) {
    return 1 / sqrt(2 * PI) * exp(-x * x / 2);
}

double *Generator1() {
    static double array[M];
    double NewSv;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.00000000001);
    for (int i = 0; i < M; i++) {
        NewSv = 0;
        for (int i = 0; i < 12; i++) {
            NewSv += dis(gen);
        }
        NewSv -= 6;
        array[i] = NewSv;
    }
    return array;
}

double *Generator2() {
    static double array[M];

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.00000000001, 1.00000000001);
    double NewGr1, NewGr2;
    for (int i = 0; i < M; i+=2) {
        NewGr1 = dis(gen); NewGr2 = dis(gen);
        array[i] = cos(2 * PI * NewGr1) * sqrt(-2 * log(NewGr2));
        array[i + 1] = sin(2 * PI * NewGr1) * sqrt(-2 * log(NewGr2));
        
    }
    return array;
}

double *Generator3() {
    static double array[M];
    double s = 0.00000000001;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.00000000001, 1.00000000001);
    double NewGr1, NewGr2;
    for (int i = 0; i < M; i += 2) {
        NewGr1 = dis(gen); NewGr2 = dis(gen);
        if (NewGr1 != NewGr2) {
            s = NewGr1 * NewGr1 + NewGr2 * NewGr2;
        }
        if (s > 1) {
            s = 1;
        }

        array[i] = NewGr1 * sqrt(-2 * log(s) / s);
        array[i + 1] = NewGr2 * sqrt(-2 * log(s) / s);
    }
    return array;
}

double* Generator4(double phi) {
    static double array[M];
    double s = 0.00000000001;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.00000000001, 1.00000000001);
    double NewGr1, NewGr2;
    for (int i = 0; i < M; i++) {
        NewGr1 = dis(gen); NewGr2 = dis(gen);
        array[i] = sqrt(-2 * log(NewGr2) * cos(2 * PI * NewGr1 - phi));
    }
    return array;
}

double* GetProb(double* L) {
    double h = (b - a) / N;
    static int array[N + 2];
    for (int i = 0; i < N + 2; i++) {
        array[i] = 0;
    }

    for (int i = 0; i < M; i++) {
        if (L[i] < a) {
            array[0]++;
        }
    }
    for (int j = 1; j < N + 1; j++) {
        for (int i = 0; i < M; i++) {
            if ((a + h * (i - 1) < L[i]) && (L[i] < a + h * i)) {
                array[j]++;
            }
        }
    }
    for (int i = 0; i < M; i++) {
        if (L[i] > b) {
            array[N + 1]++;
        }
    }

    static double PArray[N + 2];
    for (int i = 0; i < N + 2; i++) {
        PArray[i] = array[i] / M;
    }
    return PArray;
}

double GetPirs(double* P, double* v) {
    double Sum1 = 0;
    for (int i = 0; i < N + 2; i++) {
        Sum1 += (P[i] - v[i]) * (P[i] - v[i]) / v[i];
    }
    return M * Sum1;
}

int main()
{

    // Теоретические вероятности
    double ProbTheor[N+2];
    ProbTheor[0] = SympsonIntegral(-10000000, a, N);
    for (int i = 0; i < 2 * M; i++) {
        uniform_real_distribution<> dis(0.0, 1.00000000001);
    }

    // Для генератора 4
    const int k = 100;
    const double phi_0 = 0.0, phi_k = PI;

    double GenSV4[M];
    double Pirs4List[k + 1];
     

    for (int i = 0; i < k + 1; i++) {
        GetPirs(GetProb(Generator4(phi_0)),ProbTheor);
    }
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
