#include "Functions.h"

void print_arr(double** A, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cout << setw(12) << A[i][j];
        }
        cout << endl;
    }
}
void fill_random(double** A, int rows, int cols) {
    srand(time(NULL));

    int max = 10;
    int min = 1;
    double Rand = 1;

    cout << "\nEnter max for generation: ";
    enter(cin, max);
    cout << "Enter min for generation: ";
    enter(cin, min);
    cout << "Enter Rand for generation: ";
    enter(cin, Rand);

    if (min > max) {
        cout << "Error. min > max!" << endl;
        exit(991);
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            A[i][j] = (rand() % (max - min) + min) * Rand;
        }
    }
}
void fill_by_hands(double** A, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double value = 0;
            cout << "Enter A[" << i << "][" << j << "]: ";
            enter(cin, value);
            A[i][j] = value;
        }
    }
}
void sort_mattrix(double** A, int rows, int cols, int start_pos) {
    //находим строку с максимальным элементом 
    double max_elem = A[start_pos][start_pos];
    int row = start_pos;

    for (int i = start_pos; i < rows; ++i) {
        if (max_elem < A[i][start_pos]) {
            max_elem = A[i][start_pos];
            row = i;
        }
    }

    //и ставим ее в start_pos
    swap(A[row], A[start_pos]);
}
void make_triangular_matrix(double** A, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        sort_mattrix(A, rows, cols, i);//сортировка по каждому опорному элементу

        /*cout << "A_b_sort" << endl;
        print_arr(A, rows, cols);
        cout << endl;*/

        double k = A[i][i];//опорный элемент

        if (k == 0) {
            cout << "Error. Supporting element = " << k << "!" << endl;
            exit(994);
        }

        //деление строки на опорный
        for (int j = i; j < cols; ++j) {
            A[i][j] /= k;
        }

        //обнуление элементов ниже опорного
        for (int new_i = i + 1; new_i < rows; ++new_i) {
            double factor = A[new_i][i];
            for (int j = 0; j < cols; ++j) {
                A[new_i][j] -= (A[i][j] * factor);
            }
        }

        /*cout << "\nA_b_sort_null" << endl;
        print_arr(A, rows, cols);
        cout << endl;*/
    }
}

double** get_init_mattrix(int rows, int cols) {
    double** new_arr = new double* [rows];
    for (int i = 0; i < rows; ++i) {
        new_arr[i] = new double[cols] {0};
    }
    return new_arr;
}
double** get_multiplication(double** A, double** B, int A_rows, int A_cols, int B_rows, int B_cols)
{
    if (A_cols != B_rows) {
        cout << "Error. Can not to multiplicate these mattrixes!" << endl;
        exit(992);
        return nullptr;
    }

    double** result = get_init_mattrix(A_rows, B_cols);

    for (int i = 0; i < A_rows; ++i)
    {
        for (int j = 0; j < B_cols; ++j)
        {
            result[i][j] = 0;
            for (int k = 0; k < A_cols; ++k)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}
double** get_minus(double** A, double** B, int rows, int cols) {
    //должны передаваться одинаковые по размерам матрицы
    double** result = get_init_mattrix(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}
double** get_copy_mattrix(double** A, int rows, int cols) {
    double** copy_mattrix = get_init_mattrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            copy_mattrix[i][j] = A[i][j];
        }
    }
    return copy_mattrix;
}
double** get_union_mattrix(double** A, double** B, int A_rows, int A_cols, int B_rows, int B_cols) {
    if (A_rows != B_rows) {
        cout << "Error. Can not to union these mattrixs..." << endl;
        exit(993);
        return nullptr;
    }

    int un_rows = A_rows;
    int un_cols = A_cols + B_cols;
    double** un_mattrix = get_init_mattrix(un_rows, un_cols);

    //копируем А и b в un_mattrix
    for (int i = 0, i_b = 0; i < un_rows && i_b < B_rows; ++i, ++i_b) {
        for (int j = 0; j < A_cols; ++j) {
            un_mattrix[i][j] = A[i][j];
        }
        for (int j = A_cols, j_b = 0; j < un_cols && j_b < B_cols; ++j, ++j_b) {
            un_mattrix[i][j] = B[i_b][j_b];
        }
    }
    return un_mattrix;
}
double** get_roots(double** A, int rows, int cols) {
    double** roots = get_init_mattrix(rows, 1);
    int size = 0;

    roots[0][0] = A[rows - 1][cols - 1];
    size++;

    int iter = size;
    for (int i = rows - 2; i >= 0; --i) {
        double root = A[i][cols - 1];

        for (int j = cols - 2, k = 0; j >= 0 && k < size; --j, ++k) {
            root -= A[i][j] * roots[k][0];
        }
        roots[iter][0] = root;

        iter++;
        size++;
    }

    //переворот корней; в выводе х1, х2, ...
    for (int i = 0, j = size - 1; i < j; ++i, --j) {
        swap(roots[i][0], roots[j][0]);
    }
    return roots;
}

double get_max_abs_from_mattrix(double** A, int rows, int cols) {
    double max = A[0][0];
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (max < abs(A[i][j])) {
                max = abs(A[i][j]);
            }
        }
    }
    return max;
}

double** find_roots_by_Gauss(double** A, int rows, int cols) {
    //make triangular
    for (int i = 0; i < rows; ++i) {
        sort_mattrix(A, rows, cols, i);//сортировка по каждому опорному элементу

        double k = A[i][i];//опорный элемент

        //деление строки на опорный
        for (int j = i; j < cols; ++j) {
            A[i][j] /= k;
        }

        //обнуление элементов ниже опорного
        for (int j = i + 1; j < rows; ++j) {
            for (int l = 0; l < cols; ++l) {
                A[j][l] -= (A[i][l] * A[j][i]);
            }
        }
    }

    //get roots
    double** roots = get_init_mattrix(rows, 1);
    int size = 0;

    roots[0][0] = A[rows - 1][cols - 1];
    size++;

    int iter = size;
    for (int i = rows - 2; i >= 0; --i) {
        double root = A[i][cols - 1];

        for (int j = cols - 2, k = 0; j >= 0 && k < size; --j, ++k) {
            root -= A[i][j] * roots[k][0];
        }
        roots[iter][0] = root;

        iter++;
        size++;
    }

    //переворот корней; в выводе х1, х2, ...
    for (int i = 0, j = size - 1; i < j; ++i, --j) {
        swap(roots[i][0], roots[j][0]);
    }
    return roots;
}