#include <iostream>
#include "Functions.h"
#include <cmath>
#include <fstream>

using namespace std;

double fun(const double* u, double t, int n) {
    if (!n) return -u[0] * u[1] + ((t < 1e-9) ? 0.0 : (sin(t) / t));
    else return -u[1] * u[1] + (3.125 * t) / (1 + t * t);
}

void explicitMethod(const double* u, const int size) {
    double Tau;
    double Eps = 0.001, T = 1, TauMax = 0.01;
    double tk = 0;
    double* yk = new double[size];

    for (int i = 0; i < size; ++i) {
        yk[i] = u[i];
    }

    cout << "\tt" << "\t\t\tu1" << "\t\t\tu2\n";

    int iter = 0;
    do {
        double* temp = new double[size];
        //вычислить вектор f
        for (int i = 0; i < size; ++i) {
            temp[i] = fun(u, tk, i);
        }
        //определить шаг
        if (Eps / (fabs(temp[0]) + Eps / TauMax) > Eps / (fabs(temp[1]) + Eps / TauMax)) {
            Tau = Eps / (fabs(temp[1]) + Eps / TauMax);
        }
        else {
            Tau = Eps / (fabs(temp[0]) + Eps / TauMax);
        }
        //выполнить шаг
        for (int i = 0; i < size; ++i) {
            yk[i] += Tau * fabs(temp[i]);
        }
        tk += Tau;
        cout << tk << "\t\t" << yk[0] << "\t\t" << yk[1] << '\n';

        ++iter;
    } while (tk < T);

    cout << "\niterations: " << iter;

    delete[] yk;
}
void implicitMethod(const double* u, int n) {
    const double Eps = 0.001;
    double Tau, Tau_minus, Tau_plus;
    double T = 1, TauMax = 0.01, TauMin = 0.00001;
    double tk = 0, tk_plus;
    double Eps_k;
    double* yk = new double[n];
    double* yk_minus = new double[n];
    double* yk_plus = new double[n];

    Tau_minus = Tau = TauMin;

    for (int i = 0; i < n; i++) {
        yk[i] = yk_minus[i] = yk_plus[i] = u[i];
    }

    int sposob;
    cout << "\n1 - kvazioptimum, 2 - '3' zon" << endl;
    cin >> sposob;

    int kol = 0;
    do {
        do {
            tk_plus = tk + Tau;
            //решить метом ньютона систему в общем случае для tау
            yk_plus = newton(yk_plus, yk, tk, Tau, n);

            //вычислить эпсилон по формуле
            for (int k = 0; k < n; k++) {
                Eps_k = -(Tau / (Tau + Tau_minus)) * (yk_plus[k] - yk[k] - Tau * (yk[k] - yk_minus[k]) / Tau_minus);
            }
            //если модуль больше то меняем значения и переход
            for (int k = 0; k < n; k++)
                if (fabs(Eps_k) > Eps) {
                    Tau /= 2;
                    tk_plus = tk;
                    for (int j = 0; j < n; j++) {
                        yk_plus[j] = yk[j];
                    }
                }
        } while (fabs(Eps_k) > Eps);
        //определение шага согласно формулам в зависимости от выбранного варианта
        if (sposob == 1)          // квазиоптимальный
            Tau_plus = sqrt(Eps / abs(Eps_k)) * Tau;

        else if (sposob == 2)     // трехзонный
        {
            for (int i = 0; i < n; i++) {
                if (fabs(Eps_k) > Eps)
                    Tau_plus = Tau / 2;
                if ((Eps / 4 < fabs(Eps_k)) && (fabs(Eps_k) <= Eps))
                    Tau_plus = Tau;
                if (fabs(Eps_k) <= Eps / 4)
                    Tau_plus = 2 * Tau;
            }
        }
        else exit(-1);
        //если тау больше макса, то поменять местами
        if (Tau_plus > TauMax)
            Tau_plus = TauMax;
        //вывести на печать
        for (int i = 0; i < n; i++) {
            if (i == 0)
                cout << tk << setw(15);
            cout << setw(15) << yk[i];
            if (i == n - 1)
                cout << endl;
        }
        //выполнить сдвиг переменных и шагов интегрирования
        for (int i = 0; i < n; i++) {
            yk_minus[i] = yk[i];
            yk[i] = yk_plus[i];
        }
        Tau_minus = Tau;
        Tau = Tau_plus;
        tk = tk_plus;

        kol++;
    } while (tk < T);

    cout << "\niterations: " << kol << "\n";
}

int main() {
    int size = 2;
    double* y = new double[size];

    y[0] = 0.0;
    y[1] = -0.412;
    cout << "explicitMethod" << endl;
    explicitMethod(y, size);
    cout << "implicitMethod" << endl;
    implicitMethod(y, size);
}