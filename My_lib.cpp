#include"My_lib.h"

//Общие функции для матриц
double** get_init_mattrix(int rows, int cols) {
    double** new_arr = new double* [rows];
    for (int i = 0; i < rows; ++i) {
        new_arr[i] = new double[cols] {0};
    }
    return new_arr;
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


/*Метод Гаусса решения системы линейных уравнений*/
void sort_mattrix(double** A, int rows, int cols, int start_pos) {
    //находим строку с максимальным элементом 
    double max_elem = A[start_pos][start_pos];
    int row = start_pos;

    for (int i = start_pos; i < rows; ++i) {
        if (max_elem < abs(A[i][start_pos])) {
            max_elem = abs(A[i][start_pos]);
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
double** get_roots_by_Gauss(double** A, int rows, int cols) {
    //make triangular
    for (int i = 0; i < rows; ++i) {
        sort_mattrix(A, rows, cols, i);//сортировка по каждому опорному элементу

        double k = A[i][i];//опорный элемент

        //деление строки на опорный
        for (int j = i; j < cols; ++j) {
            A[i][j] /= k;
        }

        /*cout << "Del" << endl;
        print_arr(A, rows, cols);
        cout << endl;*/

        //обнуление элементов ниже опорного
        for (int j = i + 1; j < rows; ++j) {
            double d = A[j][i];
            for (int l = 0; l < cols; ++l) {
                A[j][l] -= (A[i][l] * d);
            }
        }

        /*cout << "null" << endl;
        print_arr(A, rows, cols);
        cout << endl;*/
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


/*Метод Ньютона решения системы однородных уравнений (2 уравнения с 2-мя неизвестными)*/

//Здесь по отпределению берется производная, менять не надо

double Get_dy1dx1(double x1, double x2, double M, double (*get_y1)(double, double))
{
    double y1_ = get_y1(x1 + x1 * M, x2);// + 1e-9
    double y1 = get_y1(x1, x2);

    return -1 * ((y1_ - y1) / (x1 * M));
}
double Get_dy1dx2(double x1, double x2, double M, double (*get_y1)(double, double))
{
    double y1_ = get_y1(x1, x2 + x2 * M);
    double y1 = get_y1(x1, x2);

    return -1 * ((y1_ - y1) / (x2 * M));
}
double Get_dy2dx1(double x1, double x2, double M, double (*get_y2)(double, double))
{
    double y2_ = get_y2(x1 + x1 * M, x2);
    double y2 = get_y2(x1, x2);

    return -1 * ((y2_ - y2) / (x1 * M));
}
double Get_dy2dx2(double x1, double x2, double M, double (*get_y2)(double, double))
{
    double y2_ = get_y2(x1, x2 + x2 * M);
    double y2 = get_y2(x1, x2);

    return -1 * ((y2_ - y2) / (x2 * M));
}

double** get_roots_by_Newton(double x1, double x2, double eps1, double eps2, double M, double (*get_y1)(double, double), double (*get_y2)(double, double))
{
    //Инициализация результата
    int result_rows = 2;
    int result_cols = 1;
    double** result = get_init_mattrix(result_rows, result_cols);

    //Ограничения
    const int MAX_ITER = 50;
    int k = 1;

    //инициализация Y
    int Y_rows = 2;
    int Y_cols = 1;
    double** Y = get_init_mattrix(Y_rows, Y_cols);
    Y[0][0] = get_y1(x1, x2);
    Y[1][0] = get_y2(x1, x2);

    //Инициализация счетчиков погрешностей
    double beta1 = 1;
    double beta1_2 = 1;
    double beta2_2 = 1;
    double beta2 = 1;

    //инициализация Якоби
    int J_rows = Y_rows;
    int J_cols = J_rows;
    double** J = get_init_mattrix(J_rows, J_cols);

    //вектор поправки
    double** roots = get_init_mattrix(Y_rows, Y_cols);

    while (k < MAX_ITER && beta1 >= eps1 && beta2 >= eps2) {
        //вектор невязки
        Y[0][0] = get_y1(x1, x2);
        Y[1][0] = get_y2(x1, x2);

        //Якоби
        J[0][0] = Get_dy1dx1(x1, x2, M, get_y1);
        J[0][1] = Get_dy1dx2(x1, x2, M, get_y1);
        J[1][0] = Get_dy2dx1(x1, x2, M, get_y2);
        J[1][1] = Get_dy2dx2(x1, x2, M, get_y2);

        //объединение матриц
        int J_Y_rows = J_rows;
        int J_Y_cols = J_cols + Y_cols;
        double** J_Y = get_union_mattrix(J, Y, J_rows, J_cols, Y_rows, Y_cols);

        //вектор поравки
        roots = get_roots_by_Gauss(J_Y, J_Y_rows, J_Y_cols);

        //old X
        double old_x1 = x1;
        double old_x2 = x2;

        //new X
        x1 += roots[0][0];
        x2 += roots[1][0];

        //new Y
        Y[0][0] = get_y1(x1, x2);
        Y[1][0] = get_y2(x1, x2);

        //beta1
        beta1 = get_max_abs_from_mattrix(Y, Y_rows, Y_cols);

        //beta1_2 beta2_2
        beta1_2 = abs(x1) < 1 ? abs(x1 - old_x1) : abs((x1 - old_x1) / x1);
        beta2_2 = abs(x2) < 1 ? abs(x2 - old_x2) : abs((x2 - old_x2) / x2);

        //beta2
        beta2 = beta1_2 > beta2_2 ? beta1_2 : beta2_2;

        //Result
        result[0][0] = x1;
        result[1][0] = x2;
        ++k;
    }
    if (k == MAX_ITER) {
        cout << "Error count of iteration is out of range..." << endl;
        return nullptr;
    }

    return result;
}


/*Метод прямоугольников и Симпсона на сгущающихся сетках для вычисления определенного интеграла*/
//Error - любой ссылочный параметр для записи туда погрешности
double get_integral_by_trapezoid(double (*get_func)(double), double a, double b, double eps, double& Error) {
    double n = 1;
    short p = 2;
    double h = (b - a) / n;

    double I = 1;
    double I_previous = 0;

    while (abs(I - I_previous) >= 3 * eps) {
        I_previous = I;

        double sum = 0;
        for (int i = 1; i < n; ++i) {
            sum += 2 * get_func(i * h + a);
        }

        I = ((sum + get_func(a) + get_func(b)) * h / 2);
        Error = abs((I_previous - I) / (pow(0.5, p) - 1));

        n *= 2;
        h /= 2;
    }
    //cout << "I_previous - I: " << abs(I - I_previous) << endl;
    cout << "step: " << h << endl;
    cout << "steps count: " << n << endl;

    return I;
}
double get_integral_by_Simpson(double (*get_func)(double), double a, double b, double eps, double& Error) {
    double n = 1;
    short p = 4;
    double h = (b - a) / n;

    double I = 1;
    double I_previous = 0.0;

    while (abs(I - I_previous) >= 15 * eps) {
        I_previous = I;
        double sum = 0;
        for (int i = 1; i < n; ++i) {
            sum += i % 2 == 0 ? 2 * get_func(i * h + a) : 4 * get_func(i * h + a);
        }

        I = ((sum + get_func(a) + get_func(b)) * h / 3);
        Error = abs((I_previous - I) / (pow(0.5, p) - 1));
        n *= 2;
        h /= 2;
    }
    //cout << "I_previous - I: " << abs(I - I_previous) << endl;
    cout << "step: " << h << endl;
    cout << "steps count: " << n << endl;

    return I;
}