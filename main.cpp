#include <cmath>
#include <fstream>

using namespace std;

const int NXB = 15;
const int NX = NXB * 3 + 1;
const int NYB = 12;
const int NY = NYB * 3 + 1;
const int REP = 3000;
const double EPSL = 1.e-5;
const double LL = 1.7f;
const double TEM1 = 5.0f;
const double TEM2 = 15.0f;
const double HX = 0.2f;
const double HY = 0.3f;

// расчет максимального относительного изменения температуры на текущей итерации ПВР-метода
void maxpvr(double* t1, double* del, double* maxdel) {
	double d = fabs(*del) / fabs(*t1);
	if (d > * maxdel) {
		*maxdel = d;
	}
}

int main() {
	ofstream foutT("d:/BSU/2 Cource/3 Semester/Numerical Methods/Lab_6/lab_6/dT.dat", ios_base::out | ios_base::trunc | ios_base::binary);

	double T1 = TEM1;
	double T2 = TEM2;
	double h = HX;
	double r = HY;
	double tx = 0;
	double t0 = 0;
	double t1 = 0;
	double del = 0;
	double maxdel = 0.0f;

	// содержит искомое распределение
	double** T = new double* [NY];
	for (int i0 = 0; i0 < NY; i0++) {
		T[i0] = new double[NX];
	}

	double lam = LL; // итерационный параметр λ
	double eps = EPSL;
	int prz = 1; // условие окончания ПВР-итераций
	int nT = 0; // число итераций ПВР
	// коэффициенты, используемые на итерациях ПВР
	double alf_1 = -h / r;
	double alf_2 = -r / h;
	double alf_3 = alf_2 * 0.5f;
	double alf_4 = alf_1 * 0.5f;
	double gam_1 = -2.f * (alf_1 + alf_2);
	double gam_2 = -1.5f * (alf_1 + alf_2);
	double gam_3 = -(alf_1 + alf_2);
	double gam_4 = -(alf_3 + alf_4);
	int i1 = NXB;
	int i2 = i1 + NXB;
	int i3 = i2 + NXB;
	int j1 = NYB;
	int j2 = j1 + NYB;
	int j3 = j2 + NYB;
	int rp = REP;

	// обнуление
	for (int i = 0; i <= j3; ++i) {
		for (int j = 0; j <= i3; ++j) {
			T[i][j] = 0.0f;
		}
	}

	// температуры контактов+
	for (int j = 0; j <= j1; ++j) {
		T[j][0] = T1;
	}
	for (int i = 0; i <= i1; ++i) {
		T[0][i] = T1;
	}
	for (int i = i2; i <= i3; ++i) {
		T[0][i] = T2;
	}
	for (int j = 0; j <= j1; ++j) {
		T[j][i3] = T2;
	}
	
	int k = 0;
	while (k < rp && prz == 1) {
		k++;
		for (int i = 0; i <= i3; ++i) {
			for (int j = 0; j <= j3; ++j) {
				//lines
				t0 = T[j][i];
				
				if (j == j1 && i > 0 && i < i1) { //BC+
					tx = -(alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j - 1][i]) / gam_3;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				if (i == i1 && j > j1 && j < j2) { //CD+
					tx = -(alf_4 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * T[j][i + 1]) / gam_3;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				if (j == j2 && i > i1 && i < i2) { //DE+
					tx = -(alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j - 1][i]) / gam_3;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				if (i == i2 && j > j2 && j < j3) { //EF+
					tx = -(alf_4 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * T[j][i + 1]) / gam_3;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				if (j == j3 && i > i2 && i < i3) { //FG+
					tx = -(alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j - 1][i]) / gam_3;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				if (i == i3 && j > j1 && j < j3) { //GH+
					tx = -(alf_4 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * T[j][i - 1]) / gam_3;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				if (j == 0 && i > i1 && i < i2) { //IL+
					tx = -(alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j + 1][i]) / gam_3;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				


				//points
				if (i == i1 && j == j1) {//C+
					tx = -(alf_1 * T[j - 1][i] + alf_3 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_4 * T[j + 1][i]) / gam_2;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				else if (i == i1 && j == j2) {//D+
					tx = -(alf_4 * T[j - 1][i] + alf_3 * T[j][i + 1]) / gam_4;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				else if (i == i2 && j == j2) {//E+
					tx = -(alf_1 * T[j - 1][i] + alf_3 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_4 * T[j + 1][i]) / gam_2;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				else if (i == i2 && j == j3) {//F+
					tx = -(alf_4 * T[j - 1][i] + alf_3 * T[j][i + 1]) / gam_4;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				else if (i == i3 && j == j3) {//G+
					tx = -(alf_4 * T[j - 1][i] + alf_3 * T[j][i - 1]) / gam_4;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}



				//fill
				else if (i > 0 && i < i3 && j > 0 && j < j1) {
					tx = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) / gam_1;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				else if (i > i1 && i < i3 && j > j1 - 1 && j < j2) {
					tx = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) / gam_1;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
				else if (i > i2 && i < i3 && j > j2 - 1 && j < j3) {
					tx = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) / gam_1;
					del = lam * (tx - t0); t1 = t0 + del; T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
			}
		}
		nT++;
		double w = maxdel;
		// в файл dT.dat дописывается текущая максимальная относительная поправка температуры
		foutT.write((char*)&w, sizeof w);
		if (maxdel < eps) {
			prz = 0;
		}
		maxdel = 0.0f;
	}

	foutT.close();

	// в файл nT.dat записывается число итераций ПВР
	ofstream fouT("d:/BSU/2 Cource/3 Semester/Numerical Methods/Lab_6/lab_6/nT.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	fouT.write((char*)&nT, sizeof nT);
	fouT.close();

	// в файл Pole.dat записывается распределение температуры
	ofstream fout("d:/BSU/2 Cource/3 Semester/Numerical Methods/Lab_6/lab_6/Pole.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	for (int i = 0; i < NY; ++i) {
		for (int j = 0; j < NX; ++j) {
			double w = T[i][j];
			fout.write((char*)&w, sizeof w);
		}
	}
	fout.close();

	int n_x = NX;
	int n_y = NY;

	// в файл Param.dat записывается размерность разностной сетки
	ofstream fou("d:/BSU/2 Cource/3 Semester/Numerical Methods/Lab_6/lab_6/Param.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	fou.write((char*)&n_x, sizeof n_x);
	fou.write((char*)&n_y, sizeof n_y);
	fou.close();

	for (int i0 = 0; i0 < NY; i0++) {
		delete[] T[i0];
	}
	delete[] T;

	return 0;
}
