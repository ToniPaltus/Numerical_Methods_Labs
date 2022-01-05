#pragma once
#include <iostream>
#include <istream>
#include <climits>
#include <ctime>
#include <iomanip>
#include <cmath>

using namespace std;

template <class T>
void enter(istream& in, T& data) {
    while (true) {
        in >> data;
        if (in.fail()) {
            in.clear();
            in.ignore(INT_MAX, '\n');

            cout << "An incorrect value. Please try again:";
            continue;
        }
        else {
            break;
        }
    }
}

//����� ������� ��� ������
double** get_init_mattrix(int rows, int cols);
double** get_copy_mattrix(double** A, int rows, int cols);
double** get_multiplication(double** A, double** B, int A_rows, int A_cols, int B_rows, int B_cols);
double** get_union_mattrix(double** A, double** B, int A_rows, int A_cols, int B_rows, int B_cols);
double** get_minus(double** A, double** B, int rows, int cols);

void print_arr(double** A, int rows, int cols);
void fill_random(double** A, int rows, int cols);
void fill_by_hands(double** A, int rows, int cols);

double get_max_abs_from_mattrix(double** A, int rows, int cols);


/*����� ������ ������� ������� �������� ���������*/
void sort_mattrix(double** A, int rows, int cols, int start_pos);
void make_triangular_matrix(double** A, int rows, int cols);
double** get_roots(double** A, int rows, int cols);
double** get_roots_by_Gauss(double** A, int rows, int cols);


/*����� ������� ������� ������� ���������� ��������� (2 ��������� � 2-�� ������������)*/
//����� �� ������������ ������� �����������, ������ �� ����
double Get_dy1dx1(double x1, double x2, double M, double (*get_y1)(double, double));
double Get_dy1dx2(double x1, double x2, double M, double (*get_y1)(double, double));
double Get_dy2dx1(double x1, double x2, double M, double (*get_y2)(double, double));
double Get_dy2dx2(double x1, double x2, double M, double (*get_y2)(double, double));

double** get_roots_by_Newton(double x1, double x2, double eps1, double eps2, double M, double (*get_y1)(double, double), double (*get_y2)(double, double));


/*����� ��������������� � �������� �� ����������� ������ ��� ���������� ������������� ���������*/
double get_integral_by_trapezoid(double (*get_func)(double), double a, double b, double eps, double& Error);
double get_integral_by_Simpson(double (*get_func)(double), double a, double b, double eps, double& Error);