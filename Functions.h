#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double** CreatingMatrix(int n, int m);
double** CreatingMatrix(int n);

void copyMatrix(double** X, double** copyX, int x);
void DeletingMatrix(double** X, int x);

void MatrixMultiply(double** matrix, double* roots, double* res, int size);
bool Gauss(double* An, double** X, int x);
double* newton(double* yk_plus, double* yk, double tk, double Tau, int n);

typedef double(*pf)(double*, double*, double, double);
double Differential(pf f, double* uk1, double* uk, double t, double Tau, int n);