#include"..\..\..\My_lib\My_lib\My_lib.cpp"

double get_y1(double x1, double x2) {
	return -1 * (2 * pow(x1, 2) - x1 * x2 - 5 * x1 + 1);
}
double get_y2(double x1, double x2) {
	return -1 * (x1 + 3 * log10(x1) - pow(x2, 2));
}

double get_dy1dx1(double x1, double x2) {
	return 4 * x1 - x2 - 5;
}
double get_dy1dx2(double x1, double x2) {
	return -x1;
}
double get_dy2dx1(double x1, double x2) {
	return 1 + (3 / (log(10) * x1));
}
double get_dy2dx2(double x1, double x2) {
	return -2 * x2;
}

double _Get_dy1dx1(double x1, double x2, double M) {
	double y1_ = get_y1(x1 + 1e-9, x2);// + 1e-9
	double y1 = get_y1(x1, x2);

	return -1 * ((y1_ - y1) / (1e-9));// /1e-9
}
double _Get_dy1dx2(double x1, double x2, double M) {
	double y1_ = get_y1(x1, x2 + 1e-9);
	double y1 = get_y1(x1, x2);

	return -1 * ((y1_ - y1) / (1e-9));
}
double _Get_dy2dx1(double x1, double x2, double M) {
	double y2_ = get_y2(x1 + 1e-9, x2);
	double y2 = get_y2(x1, x2);

	return -1 * ((y2_ - y2) / (1e-9));
}
double _Get_dy2dx2(double x1, double x2, double M) {
	double y2_ = get_y2(x1, x2 + 1e-9);
	double y2 = get_y2(x1, x2);

	return -1 * ((y2_ - y2) / (1e-9));
}

void print_footer(int iter, double beta1, double beta2) {
	cout << "\n\n\t\t\t|Iter: " << iter << " |Beta 1: " << setw(12) << beta1 << " |Beta 2: " << setw(12) << beta2 << endl;
}

int main() {
	//начальные приближени€
	double x1 = 3;
	double x2 = 2; //-2

	cout << "x1: " << x1 << endl;
	cout << "x2: " << x2 << endl;
	cout << "____________________________________" << endl;

	const int MAX_ITER = 50;
	int k = 1;

	double eps1 = 1e-9;
	double eps2 = 1e-9;

	/*double M = 0;
	cout << "Enter M: ";
	enter(cin, M);*/

	//инициализаци€ Y
	int Y_rows = 2;
	int Y_cols = 1;
	double** Y = get_init_mattrix(Y_rows, Y_cols);
	Y[0][0] = get_y1(x1, x2);
	Y[1][0] = get_y2(x1, x2);

	/*cout << "\nY:" << endl;
	print_arr(Y, Y_rows, Y_cols);*/

	double beta1 = 1;
	double beta1_2 = 1;
	double beta2_2 = 1;
	double beta2 = 1;

	//инициализаци€ якоби
	int J_rows = Y_rows;
	int J_cols = J_rows;
	double** J = get_init_mattrix(J_rows, J_cols);

	//вектор поправки
	double** roots = get_init_mattrix(Y_rows, Y_cols);

	while (k < MAX_ITER && beta1 >= eps1 && beta2 >= eps2) {
		//вектор нев€зки
		Y[0][0] = get_y1(x1, x2);
		Y[1][0] = get_y2(x1, x2);

		/*cout << "\nY:" << endl;
		print_arr(Y, Y_rows, Y_cols);*/

		//якоби
		J[0][0] = get_dy1dx1(x1, x2);
		J[0][1] = get_dy1dx2(x1, x2);
		J[1][0] = get_dy2dx1(x1, x2);
		J[1][1] = get_dy2dx2(x1, x2);

		/*cout << "\nJ: " << endl;
		print_arr(J, J_rows, J_cols);*/

		//объединение матриц
		int J_Y_rows = J_rows;
		int J_Y_cols = J_cols + Y_cols;
		double** J_Y = get_union_mattrix(J, Y, J_rows, J_cols, Y_rows, Y_cols);
		cout << "J_Y:" << endl;
		print_arr(J_Y, J_Y_rows, J_Y_cols);

		//вектор поравки
		roots = get_roots_by_Gauss(J_Y, J_Y_rows, J_Y_cols);
		cout << "\nroots:" << endl;
		print_arr(roots, Y_rows, Y_cols);

		//old X
		double old_x1 = x1;
		double old_x2 = x2;

		//new X
		x1 += roots[0][0];
		x2 += roots[1][0];

		cout << "\nx1: " << x1 << endl;
		cout << "x2: " << x2 << endl;

		//new Y
		Y[0][0] = get_y1(x1, x2);
		Y[1][0] = get_y2(x1, x2);

		cout << "\nY:" << endl;
		print_arr(Y, Y_rows, Y_cols);

		//beta1
		beta1 = get_max_abs_from_mattrix(Y, Y_rows, Y_cols);

		//beta1_2 beta2_2
		beta1_2 = abs(x1) < 1 ? abs(x1 - old_x1) : abs((x1 - old_x1) / x1);
		beta2_2 = abs(x2) < 1 ? abs(x2 - old_x2) : abs((x2 - old_x2) / x2);

		//beta2
		beta2 = beta1_2 > beta2_2 ? beta1_2 : beta2_2;

		print_footer(k, beta1, beta2);
		++k;
	}
	if (k == MAX_ITER) {
		cout << "Error k == MAX_ITER" << endl;
	}

	return 0;
}