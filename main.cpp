#include"../../My_lib/My_lib/My_lib.cpp"
#define K 2.0
#define A 1.0

double get_dU1(double U2, double U3) {
    return ((K - A) / A) * U2 * U3;
}
double get_dU2(double U1, double U3) {
    return ((A + K) / K) * U1 * U3;
}
double get_dU3(double U1, double U2) {
    return ((A - K) / A) * U1 * U2;
}

/*ћетод Ќьютона решени€ системы однородных уравнений (3 уравнени€ с 2-м€ неизвестными)*/

//«десь по отпределению беретс€ производна€, мен€ть не надо
double get_dy1dx2(double x2, double x3, double M, double (*get_y1)(double, double))
{
    double y1_ = get_y1(x2 + x2 * M, x3);
    double y1 = get_y1(x2, x3);

    return -1 * ((y1_ - y1) / (x2 * M));
}
double get_dy1dx3(double x2, double x3, double M, double (*get_y1)(double, double))
{
    double y1_ = get_y1(x2, x3 + x3 * M);// + 1e-9
    double y1 = get_y1(x2, x3);

    return -1 * ((y1_ - y1) / (x3 * M));
}

double get_dy2dx1(double x1, double x3, double M, double (*get_y2)(double, double))
{
    double y2_ = get_y2(x1 + x1 * M, x3);
    double y2 = get_y2(x1, x3);

    return -1 * ((y2_ - y2) / (x1 * M));
}
double get_dy2dx3(double x1, double x3, double M, double (*get_y2)(double, double))
{
    double y2_ = get_y2(x1, x3 + x3 * M);
    double y2 = get_y2(x1, x3);

    return -1 * ((y2_ - y2) / (x3 * M));
}

double get_dy3dx1(double x1, double x2, double M, double (*get_y3)(double, double))
{
    double y3_ = get_y3(x1 + x1 * M, x2);
    double y3 = get_y3(x1, x2);

    return -1 * ((y3_ - y3) / (x1 * M));
}
double get_dy3dx2(double x1, double x2, double M, double (*get_y3)(double, double))
{
    double y3_ = get_y3(x1, x2 + x2 * M);
    double y3 = get_y3(x1, x2);

    return -1 * ((y3_ - y3) / (x2 * M));
}

double** get_roots_by_Newton3(double x1, double x2, double x3, double eps1, double eps2, double M, double (*get_y1)(double, double), double (*get_y2)(double, double), double (*get_y3)(double, double))
{
    //»нициализаци€ результата
    int result_rows = 3;
    int result_cols = 1;
    double** result = get_init_mattrix(result_rows, result_cols);

    //ќграничени€
    const int MAX_ITER = 50;
    int k = 1;

    //инициализаци€ Y
    int Y_rows = 3;
    int Y_cols = 1;
    double** Y = get_init_mattrix(Y_rows, Y_cols);
    Y[0][0] = get_y1(x2, x3);
    Y[1][0] = get_y2(x1, x3);
    Y[2][0] = get_y3(x1, x2);

    //»нициализаци€ счетчиков погрешностей
    double beta1 = 1;
    double beta1_2 = 1;
    double beta2_2 = 1;
    double beta3_2 = 1;
    double beta2 = 1;

    //инициализаци€ якоби
    int J_rows = Y_rows;
    int J_cols = J_rows;
    double** J = get_init_mattrix(J_rows, J_cols);

    //вектор поправки
    double** roots = get_init_mattrix(Y_rows, Y_cols);

    while (k < MAX_ITER && beta1 >= eps1 && beta2 >= eps2) {
        //вектор нев€зки
        Y[0][0] = get_y1(x2, x3);
        Y[1][0] = get_y2(x1, x3);
        Y[2][0] = get_y3(x1, x2);

        //якоби
        J[0][0] = 0;
        J[0][1] = get_dy1dx2(x2, x3, M, get_y1);
        J[0][2] = get_dy1dx3(x2, x3, M, get_y1);

        J[1][0] = get_dy2dx1(x1, x3, M, get_y2);
        J[1][1] = 0;
        J[1][2] = get_dy2dx3(x1, x3, M, get_y2);

        J[2][0] = get_dy3dx1(x1, x2, M, get_y3);
        J[2][1] = get_dy3dx2(x1, x2, M, get_y3);
        J[2][2] = 0;

        //объединение матриц
        int J_Y_rows = J_rows;
        int J_Y_cols = J_cols + Y_cols;
        double** J_Y = get_union_mattrix(J, Y, J_rows, J_cols, Y_rows, Y_cols);

        //вектор поравки
        roots = get_roots_by_Gauss(J_Y, J_Y_rows, J_Y_cols);

        //old X
        double old_x1 = x1;
        double old_x2 = x2;
        double old_x3 = x3;

        //new X
        x1 += roots[0][0];
        x2 += roots[1][0];
        x3 += roots[2][0];

        //new Y
        Y[0][0] = get_y1(x2, x3);
        Y[1][0] = get_y2(x1, x3);
        Y[2][0] = get_y2(x1, x2);

        //beta1
        beta1 = get_max_abs_from_mattrix(Y, Y_rows, Y_cols);
         
        //beta1_2 beta2_2 beta3_2
        beta1_2 = abs(x1) < 1 ? abs(x1 - old_x1) : abs((x1 - old_x1) / x1);
        beta2_2 = abs(x2) < 1 ? abs(x2 - old_x2) : abs((x2 - old_x2) / x2);
        beta3_2 = abs(x3) < 1 ? abs(x3 - old_x3) : abs((x3 - old_x3) / x3);

        //beta2
        double max = beta1_2;
        if (max < beta2_2) {
            max = beta2_2;
        }
        if (max < beta3_2) {
            max = beta3_2;
        }

        beta2 = max;

        //Result
        result[0][0] = x1;
        result[1][0] = x2;
        result[2][0] = x3;
        ++k;
    }
    if (k == MAX_ITER) {
        cout << "Error count of iteration is out of range..." << endl;
        return nullptr;
    }

    return result;
}

double get_min_from(double a, double b, double c) {
    double min = a;
    if (min > b) {
        min = b;
    }
    if (min > c) {
        min = c;
    }
    return min;
}
double get_max_from(double a, double b, double c) {
    double max = a;
    if (max < b) {
        max = b;
    }
    if (max < c) {
        max = c;
    }
    return max;
}

//3 вариант
double** get_solution_by_explicit_method(double t, double T, double EPS, double t_max, double U1, double U2, double U3) {
	double** result = get_init_mattrix(3, 1);
	double h = 0;
	int i = 1;
	while (t < T) {
        double dU1 = get_dU1(U2, U3);
		double dU2 = get_dU2(U1, U3);
        double dU3 = get_dU3(U1, U2);

		//max U1 U2 U3 = max
		double max = U1;
		if (U2 > max) {
			max = U2;
		}
		if (U3 > max) {
			max = U3;
		}

		double eps = EPS * abs(max);

		double h1 = eps / (abs(dU1) + eps / t_max);
		double h2 = eps / (abs(dU2) + eps / t_max);
		double h3 = eps / (abs(dU3) + eps / t_max);

		//find min h1 h2 h3 = h
		h = h1;
		if (h2 < h) {
			h = h2;
		}
		if (h3 < h) {
			h = h3;
		}

		//update U1 U2 U3
		U1 += h * dU1;
		U2 += h * dU2;
		U3 += h * dU3;

		t += h;
		cout << "Iter: " << i++ << endl;

		//data
		cout << setw(6) << "U1: " << setw(6) << U1 << ";" << setw(6)
			<< setw(6) << "U2: " << setw(6) << U2 << ";" << setw(6)
			<< setw(6) << "U3: " << setw(6) << U3 << ";" << setw(6)
			<< setw(6) << "h: " << setw(6) << h << ";" << setw(6)
            << setw(6) << "t: " << setw(6) << t << ";\n" << endl;
	}
	result[0][0] = U1;
	result[1][0] = U2;
	result[2][0] = U3;

	return result;
}


int main() {
    /*
	double** solutions_explicit = get_solution_by_explicit_method(0, 1, 1e-2, 1e-3, 1, 1, 1);
	cout << "Explicit method:" << endl;
	print_arr(solutions_explicit, 3, 1);
    */

    double U1 = 1;
    double U2 = 1;
    double U3 = 1;

    double T = 1;
    double EPS = 1e-3;
    double tau_min = 1e-3;
    double tau_max = 1e-2;
    
    double t = 0;
    double t_next = 0;

    double U1_prev, U1_cur, U1_next = 0;
    double U2_prev, U2_cur, U2_next = 0;
    double U3_prev, U3_cur, U3_next = 0;

    U1_prev = U1_cur = U1_next = U1;
    U2_prev = U2_cur = U2_next = U2;
    U3_prev = U3_cur = U3_next = U3;

    double tau_prev, tau_cur = 0;
    tau_prev = tau_cur = tau_min;

    int i = 0;
    while (t < T) {
        t_next = t + tau_cur;

        double** new_dU = get_roots_by_Newton3(U1_cur, U2_cur, U3_cur, 1e-9, 1e-9, 1e-9, get_dU1, get_dU2, get_dU3);
        //print_arr(new_U, 3, 1);
        U1_next = new_dU[0][0];
        U2_next = new_dU[1][0];
        U3_next = new_dU[2][0];

        // EPS = 0.01 * abs(get_max_from(U1_next, U2_next, U3_next));

        double eps1 = -(tau_cur/(tau_cur + tau_prev))*(U1_next - U1_cur - (tau_cur / tau_prev) * (U1_cur - U1_prev));
        double eps2 = -(tau_cur/(tau_cur + tau_prev))*(U2_next - U2_cur - (tau_cur / tau_prev) * (U2_cur - U2_prev));
        double eps3 = -(tau_cur/(tau_cur + tau_prev))*(U3_next - U3_cur - (tau_cur / tau_prev) * (U3_cur - U3_prev));

        if (abs(eps1) > EPS) {
            tau_cur /= 2;
            t_next = t;
            U1_next = U1_cur;
            continue;
        }
        if (abs(eps2) > EPS) {
            tau_cur /= 2;
            t_next = t;
            U2_next = U2_cur;
            continue;
        }
        if (abs(eps3) > EPS) {
            tau_cur /= 2;
            t_next = t;
            U3_next = U3_cur;
            continue;
        }


        double tau1_next = 0;
        if (abs(eps1) > EPS) {
            tau1_next = tau_cur / 2;
        }
        if ((abs(eps1) > EPS / 4) && (abs(eps1) <= EPS)) {
            tau1_next = tau_cur;
        }
        if (abs(eps1) <= EPS / 4) {
            tau1_next = 2 * tau_cur;
        }

        double tau2_next = 0;
        if (abs(eps2) > EPS) {
            tau2_next = tau_cur / 2;
        }
        if ((abs(eps2) > EPS / 4) && (abs(eps2) <= EPS)) {
            tau2_next = tau_cur;
        }
        if (abs(eps2) <= EPS / 4) {
            tau2_next = 2 * tau_cur;
        }

        double tau3_next = 0;
        if (abs(eps3) > EPS) {
            tau3_next = tau_cur / 2;
        }
        if ((abs(eps3) > EPS / 4) && (abs(eps3) <= EPS)) {
            tau3_next = tau_cur;
        }
        if (abs(eps3) <= EPS / 4) {
            tau3_next = 2 * tau_cur;
        }

        //очень маленький тау
        double tau_next = get_min_from(tau1_next, tau2_next, tau3_next);
        if (tau_next > tau_max) {
            t_next = tau_max;
        }
        if (t_next < tau_min) {
            t_next = tau_min;
        }

        cout << "Iter: " << i++ << endl <<
            "U1_next: " << U1_next << "| U2_next: " << U2_next <<
            "| U3_next: " << U3_next << "| t_next: " << t_next << endl << endl;

        U1_prev = U1_cur;
        U2_prev = U2_cur;
        U3_prev = U3_cur;

        U1_cur = U1_next;
        U2_cur = U2_next;
        U3_cur = U3_next;

        tau_prev = tau_cur;
        tau_cur = tau_next;
        t = t_next;
    }
    
	return  0;
}