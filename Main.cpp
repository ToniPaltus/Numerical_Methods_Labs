/* Задание
 Решить систему линейных уравнений методом Гаусса и посчитать погрешности копьютерных вычислений
 */
#include "Functions.h"

int main() {
    //A
    int A_rows = 0;
    cout << "Enter A_rows: ";
    enter(cin, A_rows);
    int A_cols = A_rows;

    double** A = get_init_mattrix(A_rows, A_cols);
    //fill_random(A, A_rows, A_cols);
    fill_by_hands(A, A_rows, A_cols);
    cout << "\nA:\n" << endl;
    print_arr(A, A_rows, A_cols);

    //b
    int b_rows = A_rows;
    int b_cols = 1;

    double** b = get_init_mattrix(b_rows, b_cols);
    //fill_random(b, b_rows, b_cols);
    fill_by_hands(b, b_rows, b_cols);
    cout << "\nb:\n" << endl;
    print_arr(b, b_rows, b_cols);

    //union A and b
    int A_b_rows = A_rows;
    int A_b_cols = A_cols + b_cols;

    double** A_b = get_union_mattrix(A, b, A_rows, A_cols, b_rows, b_cols);
    cout << "\nA_b:\n" << endl;
    print_arr(A_b, A_b_rows, A_b_cols);

    //gauss->
    //triangular A_b
    make_triangular_matrix(A_b, A_b_rows, A_b_cols);
    cout << "\nA_b triangular:\n" << endl;
    print_arr(A_b, A_b_rows, A_b_cols);

    //roots A_b
    int A_b_roots_rows = A_b_rows;
    int A_b_roots_cols = 1;

    double** A_b_roots = get_roots(A_b, A_b_rows, A_b_cols);
    cout << "\nA_b_roots:\n" << endl;
    print_arr(A_b_roots, A_b_roots_rows, A_b_roots_cols);

    cout << "\ntestGausstestGausstestGausstestGausstestGausstestGausstestGauss";
    double** test = find_roots_by_Gauss(A_b, A_b_rows, A_b_cols);
    cout << "\ntest:" << endl;
    print_arr(test, A_b_rows, 1);
    cout << "testGausstestGausstestGausstestGausstestGausstestGausstestGauss" << endl;

    //tests ->
    //new b
    int new_b_rows = A_rows;
    int new_b_cols = A_b_roots_cols;

    double** new_b = get_multiplication(A, A_b_roots, A_rows, A_cols, A_b_roots_rows, A_b_roots_cols);
    cout << "\nnew_b:\n" << endl;
    print_arr(new_b, new_b_rows, new_b_cols);

    //minus F_1
    int F_1_rows = new_b_rows;
    int F_1_cols = new_b_cols;
    double** F_1 = get_minus(new_b, b, F_1_rows, F_1_cols);
    cout << "\nnew_b - b:\n" << endl;
    print_arr(F_1, F_1_rows, F_1_cols);

    //max from F_1
    double delta_1 = get_max_abs_from_mattrix(F_1, F_1_rows, F_1_cols);
    cout << "\ndelta_1: " << delta_1 << endl;

    //union A and new_b
    int A_new_b_rows = A_rows;
    int A_new_b_cols = A_cols + new_b_cols;
    double** A_new_b = get_union_mattrix(A, new_b, A_rows, A_cols, new_b_rows, new_b_cols);
    cout << "\nA_new_b:\n" << endl;
    print_arr(A_new_b, A_new_b_rows, A_new_b_cols);

    //triangular A_new_b
    make_triangular_matrix(A_new_b, A_new_b_rows, A_new_b_cols);
    cout << "\nA_new_b triangular:\n" << endl;
    print_arr(A_new_b, A_new_b_rows, A_new_b_cols);

    //roots A_new_b
    int A_new_b_roots_rows = A_new_b_rows;
    int A_new_b_roots_cols = 1;

    double** A_new_b_roots = get_roots(A_new_b, A_new_b_rows, A_new_b_cols);
    cout << "\nA_new_b_roots:\n" << endl;
    print_arr(A_new_b_roots, A_new_b_roots_rows, A_new_b_roots_cols);

    //minus F_2
    int F_2_rows = new_b_rows;
    int F_2_cols = new_b_cols;
    double** F_2 = get_minus(A_new_b_roots, A_b_roots, A_new_b_roots_rows, A_new_b_roots_cols);
    cout << "\nA_new_b_roots - A_b_roots:\n" << endl;
    print_arr(F_2, F_2_rows, F_2_cols);

    //max from F_2
    double delta_2 = get_max_abs_from_mattrix(F_2, F_2_rows, F_2_cols);
    cout << "\nmax abs from vector (A_new_b_roots - A_b_roots): " << delta_2 << endl;

    //max from A_b_roots
    double max_A_b_roots = get_max_abs_from_mattrix(A_b_roots, A_b_roots_rows, A_b_roots_cols);
    cout << "\nmax_A_b_roots: " << max_A_b_roots << endl;

    //sigma
    double sigma = delta_2 / max_A_b_roots;
    cout << "\nsigma: " << sigma << endl;

    return 0;
}