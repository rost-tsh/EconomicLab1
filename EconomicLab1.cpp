// EconomicLab1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#define _USE_MATH_DEFINES

#include <cmath>
#include <random>


using namespace std;


const double PI = 3.141592653589793;
const double A = -3, B = 3;
const int N = 1000;
const int M = 10000;
const double E = 0.00001;

double f(double x)
{
    return 1 / sqrt(2 * PI) * exp(-x * x / 2);
}

double CentralIntegral(double a, double b, double h)
{
    double Integral = 0.0;
    double n = (b - a) / h;
    for (int i = 1; i <= n; i++)
        Integral += h * f(a + h * (i - 0.5));
    return Integral;
}

double TrapezeIntegral(double a, double b, double h)
{
    double n = (b - a) / h;
    double Integral = h * (f(a) + f(b)) / 2.0;
    for (int i = 1; i <= n - 1; i++)
        Integral = Integral + h * f(a + h * i);
    return Integral;
}

double SympsonIntegral(double a, double b, double h)
{
    double n = (b - a) / h;
    double Integral = h * (f(a) + f(b)) / 6.0;
    for (int i = 1; i <= n; i++)
        Integral = Integral + 4.0 / 6.0 * h * f(a + h * (i - 0.5));
    for (int i = 1; i <= n - 1; i++)
        Integral = Integral + 2.0 / 6.0 * h * f(a + h * i);

    return Integral;
}

double p(double x) {
    return 1 / sqrt(2 * PI) * exp(-x * x / 2);
}

double MinElem(double* L, int N) {
    int min = 0;
    for (int i = 1; i < N; i++)
    {
        if (L[min] > L[i]) min = i;
    }
    return L[min];
}

double *Generator1() {
    static double array[M];
    double NewSv;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
    for (int i = 0; i < size(array); i++) {
        NewSv = 0;
        for (int j = 0; j < 12; j++) {
            NewSv += dis(gen);
        }
        NewSv -= 6;
        array[i] = NewSv;
    }
    return array;
}

double* Generator2() {
    static double array[M];

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(E, 1);
    double NewGr1, NewGr2;
    for (int i = 0; i < size(array); i+=2) {
        NewGr1 = dis(gen); NewGr2 = dis(gen);
        array[i] = cos(2 * PI * NewGr1) * sqrt(-2 * log(NewGr2));
        array[i + 1] = sin(2 * PI * NewGr1) * sqrt(-2 * log(NewGr2));
        
    }
    return array;
}

double* Generator3() {
    static double array[M];
    double s;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-1, 1);
    double NewGr1, NewGr2;
    for (int i = 0; i < size(array); i += 2) {
        s = E;
        NewGr1 = dis(gen); NewGr2 = dis(gen);
        if (abs(NewGr1 - NewGr2) >= E) {
            s = NewGr1 * NewGr1 + NewGr2 * NewGr2;
            
        }
        if (abs(s - 1) < E) {
            s = 1;
        }
        array[i] = NewGr1 * sqrt(-2 * log(s) / s);
        array[i + 1] = NewGr2 * sqrt(-2 * log(s) / s);
    }
    return array;
}

double* Generator4(double phi) {
    static double array[M];
    double s = E;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(E, 1-E);
    double NewGr1, NewGr2;
    for (int i = 0; i < size(array); i++) {
        NewGr1 = dis(gen); NewGr2 = dis(gen);
        array[i] = sqrt(-2 * log(NewGr2)) * cos(2 * PI * NewGr1 - phi);
        
    }
    return array;
}

int* GetM(double* L) {
    double h = (B - A) / N;
    static int array[N + 2];
    for (int i = 0; i < size(array); i++) {
        array[i] = 0;
    }

    for (int i = 0; i < M; i++) {
        if (L[i] <= A) {
            array[0]++;
        }
    }
    for (int j = 1; j < size(array) - 1; j++) {
        for (int i = 0; i < M; i++) {
            if ((A + h * (j - 1) <= L[i]) && (L[i] < A + h * j)) {
                array[j]++;
            }
        }
    }
    for (int i = 0; i < M; i++) {
        if (L[i] >= B) {
            array[size(array) - 1]++;
        }
    }
    return array;
}

double* GetProb(double* L) {
    int* array = GetM(L);
    static double PArray[N + 2];
    for (int i = 0; i < N + 2; i++) {
        PArray[i] = (double)array[i] / (double)(M);
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

void PrintArray(double* P, int n) {
    for (int i = 0; i < n; i++) {
        cout << P[i] << endl;
    }
}


int main()
{
    setlocale(LC_ALL, "Russian");
    double h = (B - A) / N;
    static double array[M];

    // Теоретические вероятности
    static double ProbTheor[N+2];

    ProbTheor[0] = SympsonIntegral(-100, A, 0.005);

    for (int i = 1; i < N + 1; i++) {
        ProbTheor[i] = SympsonIntegral(A + h * (i - 1), A + h * i, 0.005);
    }
    ProbTheor[N+1] = SympsonIntegral(B, 100, 0.005);

    // Для генератора 1
    double Pirs1 = GetPirs(GetProb(Generator1()), ProbTheor);
    cout << "Для генератора 1 критерий согласия Пирсона = " << Pirs1 << endl;

    // Для генератора 2
    double Pirs2 = GetPirs(GetProb(Generator2()), ProbTheor);
    cout << "Для генератора 2 критерий согласия Пирсона = " << Pirs2 << endl;

    // Для генератора 3
    double Pirs3 = GetPirs(GetProb(Generator3()), ProbTheor);
    cout << "Для генератора 3 критерий согласия Пирсона = " << Pirs3 << endl;

    // Для генератора 4
    const int k = 100;
    const double phi_0 = 0.0, phi_k = PI;
    double phi_h = (phi_k - phi_0) / k;
    double Pirs4List[k + 1];
    for (int i = 0; i < k + 1; i++) {
        Pirs4List[i] = GetPirs(GetProb(Generator4(phi_0 + phi_h * i)), ProbTheor);
    }
    double Pirs4 = MinElem(Pirs4List, k + 1);
    cout << "Для генератора 4 критерий согласия Пирсона = " << Pirs4 << endl;

   
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
