#include"..\..\My_lib\My_lib\My_lib.cpp"

double get_sum(double** A, int rows, int cols, int start_degree) {
	double sum = 0;
	for (int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j){
			sum += pow(A[i][j], start_degree);
		}
	}
	return sum;
}
double get_sum_B(double** A, double** B, int rows, int cols, int start_degree) {
	// A - D, B - V
	double sum = 0;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			sum += (B[i][j] * pow(A[i][j], start_degree));
		}
	}
	return sum;
}
void fill_A(double** A, double** D, int A_rows, int A_cols, int D_rows, int D_cols) {
	for (int i = 0; i < A_rows; ++i) {
		int degree = i;
		for (int j = 0; j < A_cols; ++j) {
			A[i][j] = get_sum(D, D_rows, D_cols, degree++);
		}
	}
}
void fill_B(double** B, double** D, double** V, int B_rows,int B_cols,int D_rows,int D_cols) {
	for (int i = 0; i < B_rows; ++i) {
		int degree = i;
		for (int j = 0; j < B_cols; ++j) {
			B[i][j] = get_sum_B(D, V, D_rows, D_cols, degree);
		}
	}
}
double get_dispersia(double** D, double** V, double** a,int cols) {
	double disp = 0;
	for (int i = 0; i < cols; ++i) {
		disp += pow(V[0][i] - a[0][0] - a[1][0] * D[0][i] - a[2][0] * pow(D[0][i], 2), 2);
	}
	disp /= (10 - 2 - 1);
	return disp;
}

double get_f(double** a, double x){
	double result = a[0][0] + a[1][0] * x + a[2][0] * pow(x, 2);
	return result;
}
double** get_points_f(double** D,double** a, int D_rows, int D_cols){
	double** result = get_init_mattrix(D_rows, D_cols);
	for (int i = 0; i < D_cols; ++i) {
		result[0][i] = get_f(a, D[0][i]);
	}
	return result;
}

int main() {
	//Начальные данные
	int D_rows = 1;
	int D_cols = 10;
	double** D = get_init_mattrix(D_rows, D_cols);

	for (int i = 0; i < D_cols; ++i) {
		D[0][i] = static_cast<double>(i) / 10;
	}
	cout << "X:" << endl;
	print_arr(D, D_rows, D_cols);

	double** V = get_init_mattrix(D_rows, D_cols);
	V[0][0] = 0.957;
	V[0][1] = 0.969;
	V[0][2] = 0.976;
	V[0][3] = 0.978;
	V[0][4] = 0.975;
	V[0][5] = 0.968;
	V[0][6] = 0.954;
	V[0][7] = 0.939;
	V[0][8] = 0.918;
	V[0][9] = 0.894;

	cout << "Y:" << endl;
	print_arr(V, D_rows, D_cols);

	//A
	int A_rows = 3;
	int A_cols = A_rows;
	double** A = get_init_mattrix(A_rows, A_cols);
	fill_A(A, D, A_rows, A_cols, D_rows, D_cols);
	/*cout << "A:" << endl;
	print_arr(A, A_rows, A_cols);*/

	//B
	int B_rows = A_rows;
	int B_cols = 1;
	double** B = get_init_mattrix(B_rows, B_cols);
	fill_B(B, D, V, B_rows, B_cols, D_rows, D_cols);
	/*cout << "B:" << endl;
	print_arr(B, B_rows, B_cols);*/

	//A_B
	int A_B_rows = A_rows;
	int A_B_cols = A_cols + B_cols;
	double** A_B = get_union_mattrix(A, B, A_rows, A_cols, B_rows, B_cols);
	cout << "A_B:" << endl;
	print_arr(A_B, A_B_rows, A_B_cols);

	//a
	int a_rows = A_rows;
	int a_cols = 1;
	double** a = get_roots_by_Gauss(A_B, A_B_rows, A_B_cols);
	cout << "a:" << endl;
	print_arr(a, a_rows, a_cols);

	//disp + sigma(среднеквадратичное отклонение)
	double disp = get_dispersia(D, V, a, D_cols);
	cout << "\nDispersia^2: " << disp << endl;
	double sigma = sqrt(disp);
	cout << "Sigma: " << sigma << endl;

	//new Y
	double** new_f = get_points_f(D, a, D_rows, D_cols);
	cout << "\nNew Y:" << endl;
	print_arr(new_f, D_rows, D_cols);

	return 0;
}