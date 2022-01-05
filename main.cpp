#include "..\..\My_lib\My_lib\My_lib.cpp"

double get_func(double x) {
	return 1 / (pow(log(x), 2) + 1);
}

//1.327865229389399
int main() {
	double a = 1;
	double b = 2.835;

	double eps = 1e-5;
	double Error = 0.0;

	double I1 = get_integral_by_trapezoid(get_func, a, b, eps, Error);
	cout << "I1 = " << I1 << " +- " << Error << endl << endl;
	double I2 = get_integral_by_Simpson(get_func, a, b, eps, Error);
	cout << "I2 = " << I2 << " +- " << Error << endl;

	return 0;
}