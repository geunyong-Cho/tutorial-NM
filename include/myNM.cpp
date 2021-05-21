/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [Cho Geunyong]
Created          : 26-03-2018
Modified         : 30-04-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"

// Original function
double original_func(double _x) {
	double L = 4; // [m]
	double E = 70e+6; // [kN / m^4]
	double I = 52.9e-6; // [m^4]
	double W_0 = 20; // [kN / m]

	double y = ((W_0 * _x) * (7 * L * L * L * L - 10 * L * L * _x * _x + 3 * _x * _x * _x * _x)) / (360 * E * I * L);
	return y;
}

// Function Differentiation
double ass1_func(double _x) {

	double y = (125000 * _x * _x * _x * _x) / 11109 - (20000000 * _x * _x) / 33327 - (125000 * _x * (-12 * _x * _x * _x + 320 * _x)) / 33327 + 32000000 / 4761;// 1 diff function y

	return y;
}

double ass1_dfunc(double _x) {
	double y = (125000 * _x * (36 * _x * _x - 320)) / 33327 - (80000000 * _x * _x) / 33327 + (1000000 * _x * _x * _x) / 11109; // 2 diff function y
	return y;
}

// Bisection method
double bisection(double _a0, double _b0, double _tol)
{
	int k = 0;
	int Nmax = 1000;
	double a = _a0;
	double b = _b0;
	double xn = 0;
	double ep = 1000;

	do {
		xn = (a + b) / 2;
		ep = fabs(ass1_func(xn));
		k++;
		printf("Iteration:%d \t", k);
		printf("X(n): %f \t", xn);
		printf("Tolerance: %.10f\n", ep);

		if (ass1_func(a) * ass1_func(xn) < 0)   // f(a)*f(xn)이 0보다 큰지 작은지에 따라 boundary 설정 다르게 함
			b = xn;
		else
			a = xn;
	} while (k<Nmax && ep>_tol);                // do-while 반복 조건
	return xn;
}

// Newton Raphson method
double newtonRaphson(double _x0, double _tol)
{
	int k = 0;
	int Nmax = 100;
	double xk;
	double xk1;
	xk = _x0;
	double ep = 0;
	double hk;
	do {
		if (ass1_dfunc(xk) != 0) {                        // f'(xk) = 0 이 되면 Newton Raphson Method 사용 불가
			hk = -1 * ass1_func(xk) / ass1_dfunc(xk);    // hk 설정, xk = xk + hk
			xk1 = xk + hk;
			ep = fabs((xk1 - xk) / xk);
			k++;
			printf("Iteration:%d \t", k);
			printf("X(k): %f \t", xk);
			printf("Tolerance: %.10f\n", ep);
			xk = xk1;
		}
		else {
			printf("Wrong starting value");
			break;
		}

	} while (k<Nmax && ep>_tol);
	return xk;
}

double secant(double _x1, double _x2, double _tol)
{
	int k = 0;
	int Nmax = 100;
	double xk = _x1;
	double xk1 = _x2;
	double xk2;
	double ep = 0;
	double dfunc;
	double hk;
	do {
		if (ass1_dfunc(xk) != 0) {
			dfunc = (ass1_func(xk1) - ass1_func(xk)) / (xk1 - xk);
			hk = -1 * ass1_func(xk1) / dfunc;
			xk2 = xk1 + hk;
			ep = fabs((xk2 - xk1) / xk1);
			k++;
			printf("Iteration:%d \t", k);
			printf("X(k): %f \t", xk1);
			printf("Tolerance: %.10f\n", ep);
			xk = xk1;
			xk1 = xk2;

		}
		else {
			printf("Wrong starting value");
			break;
		}

	} while (k<Nmax && ep>_tol);
	return xk1;
}

// Bonus Original function
double bonus_func(double _x) {
	double y = (1 / _x) - 2;
	return y;
}

// Bonus 1diff function
double bonus_dfunc(double _x) {
	double y = -1 / (_x * _x); // differentiation of y
	return y;
}
// Bonus Newton Raphson
double bonus_newtonRaphson(double _x0, double _tol)
{
	int k = 0;
	int Nmax = 100;
	double xk = _x0;
	double xk1;
	double ep = 0;
	double hk;
	do {
		if (bonus_dfunc(xk) != 0) {
			hk = -1 * bonus_func(xk) / bonus_dfunc(xk);
			xk1 = xk + hk;
			ep = fabs((xk1 - xk) / xk);
			xk = xk1;

			k++;
			printf("Iteration:%d \t", k);
			printf("X(k): %f \t", xk);
			printf("Tolerance: %.5f\n\n", ep);
		}
		else {
			printf("Wrong starting value");
			break;
		}

	} while (k<Nmax && ep>_tol);
	return xk;
}
// Hybrid method
double hybrid_method(double _a0, double _b0, double _x0, double _tol) {
	int k = 1;
	int Nmax = 100;
	double xk = _x0;
	double xk1;
	double ep = 1000;
	double hk;
	do {
		if (bonus_func(xk) != 0) {
			if ((xk > _b0) || (xk < _a0) || ep < 0.1) {          // xk가 boundary를 벗어나거나 ep가 작으면 bisection method 사용
				xk = (_a0 + _b0) / 2;
				ep = fabs(bonus_func(xk));
				if (bonus_func(_a0) * bonus_func(xk) < 0)
					_b0 = xk;
				else
					_a0 = xk;
				printf("Bisection \t");
				printf("Iteration:%d \t", k);
				printf("X(n): %f \t\t", xk);
				printf("Tolerance: %.10f\n", ep);
				k++;
			}
			else {                                               // 나머지 경우 Newton Raphson Method 사용
				hk = -1 * bonus_func(xk) / bonus_dfunc(xk);
				xk1 = xk + hk;
				ep = fabs((xk1 - xk) / xk);
				xk = xk1;
				printf("Newton Raphson \t");
				k++;
				printf("Iteration:%d \t", k);
				printf("X(k): %f \t", xk);
				printf("Tolerance: %.10f\n", ep);
			}
		}
		else {
			printf("Wrong starting value");
			break;
		}
	} while (k<Nmax && ep>_tol);
	return xk;
}

//function / dfunction
double func1(double _x, double _y)
{
	return _y - (exp(_x / 2) + exp(-_x / 2)) / 2;
}

double func2(double _x, double _y)
{
	return 9 * pow(_x, 2) + 25 * pow(_y, 2) - 225;
}

double dfunc1x(double _x, double _y)
{
	return -(exp(_x / 2) - exp(-_x / 2)) / 4;
}

double dfunc1y(double _x, double _y)
{
	return 1;
}

double dfunc2x(double _x, double _y)
{
	return 18 * _x;
}

double dfunc2y(double _x, double _y)
{
	return 50 * _y;
}

Matrix FX(double _x, double _y)
{
	Matrix F = createMat(2, 1);
	F.at[0][0] = func1(_x, _y);
	F.at[1][0] = func2(_x, _y);
	return F;
}

Matrix Jac(double _x, double _y)
{
	Matrix J = createMat(2, 2);
	J.at[0][0] = dfunc1x(_x, _y);
	J.at[0][1] = dfunc1y(_x, _y);
	J.at[1][0] = dfunc2x(_x, _y);
	J.at[1][1] = dfunc2y(_x, _y);
	return J;
}

void System_non_linear(double _x, double _y, double _tol)
{
	double xk = _x;
	double yk = _y;
	int k = 0;
	double ep = 0;
	int Nmax = 20;
	Matrix Xk = createMat(2, 1);
	Xk.at[0][0] = _x;
	Xk.at[1][0] = _y;
	Matrix Xk1 = createMat(2, 1);
	Matrix hk = createMat(2, 1);
	Matrix U = zeros(2, 2);
	Matrix d = zeros(2, 1);
	Matrix F = zeros(2, 1);
	Matrix J = zeros(2, 2);
	do
	{
		k++;
		printf("Iteration:%d \t", k);
		F = FX(xk, yk);
		J = Jac(xk, yk);
		U = gaussElimU(J, F, U, d);
		d = gaussElimd(J, F, U, d);
		/*printMat(F, "\nF");
		printMat(J, "\nJ");
		printMat(U, "\nU");
		printMat(d, "\nd");*/
		hk = backSub(U, minus(d));
		Xk1 = addMat(Xk, hk);
		/*printMat(hk, "\nhk");
		printMat(Xk, "\nXk");
		printMat(Xk1, "\nXk1");*/
		ep = fabs((normMat(Xk1) - normMat(Xk)) / normMat(Xk));
		Xk = copyMat(Xk1);
		xk = Xk.at[0][0];
		yk = Xk.at[1][0];
	} while (k < Nmax && ep > _tol);
}

// data set x, y를 받아 1차 curve fitting 실시하여 1행을 계수 a1, 2행을 계수 a0으로 갖는 행렬 z 반환
Matrix	linearFit(Matrix _x, Matrix _y) {
	int mx = _x.rows;
	int my = _y.rows;
	int m;
	m = mx;
	double a1 = 0;
	double a0 = 0;
	Matrix z = createMat(2, 1);
	if (mx != my) {
		printf("Error: Dataset should have same point.\n");
	}
	if ((mx < 2) || (my < 2)) {
		printf("Error: Dataset should be more than 2.\n");
	}
	else {
		double Sx = 0; double Sxx = 0; double Sxy = 0; double Sy = 0;
		a1 = 0; a0 = 0;
		for (int k = 0; k < m; k++) {
			Sx += _x.at[k][0];
			Sxx += _x.at[k][0] * _x.at[k][0];
			Sy += _y.at[k][0];
			Sxy += _x.at[k][0] * _y.at[k][0];
		}
		if ((m * Sxx - Sx * Sx != 0) && (m * Sxx - Sx * Sx != 0)) {
			a1 = (m * Sxy - Sx * Sy) / (m * Sxx - Sx * Sx);
			a0 = (-1 * Sx * Sxy + Sxx * Sy) / (m * Sxx - Sx * Sx);
		}
		else {
			printf("Error : 0 Division\n");
		}
	}
	z.at[0][0] = a1;
	z.at[1][0] = a0;
	printMat(_x, "T");
	printMat(_y, "P");
	printMat(z, "z");
	return z;
}

// data set x, y, interpolation을 통해 찾고자 하는 지점 xq를 입력받아 xq에 대응하는 data yq를 반환
double linearinterp(Matrix _x, Matrix _y, int xq)
{	
	int mx = _x.rows;
	int my = _y.rows;
	int index_end = _x.rows + 10;
	int intv = 5;
	double a1; double a0;
	Matrix x2 = createMat(index_end, 1);
	Matrix y2 = copyMat(x2);
	if (mx != my) {
		printf("Error: Dataset should have same point.\n");
	}
	if ((mx < 2) || (my < 2)) {
		printf("Error: Dataset should be more than 2.\n");
	}
	else {
		for (int k = 0; k < index_end; k++) {
			x2.at[k][0] = k * intv;
		}

		for (int i = 0; i < _x.rows - 1; i++) {
			if (_x.at[i + 1][0] - _x.at[i][0] != 0) {
				a1 = (_y.at[i + 1][0] - _y.at[i][0]) / (_x.at[i + 1][0] - _x.at[i][0]);
				a0 = _y.at[i][0];
				y2.at[2 * i + 1][0] = a1 * (x2.at[2 * i + 1][0] - _x.at[i][0]) + a0;
				y2.at[2 * i][0] = _y.at[i][0];
			}
			else {
				printf("Error : 0 Division\n");
			}
		}
	}
	y2.at[index_end - 1][0] = _y.at[_y.rows - 1][0];
	printMat(x2, "X2");
	printMat(y2, "Y2");

	int a = (x2.at[index_end - 1][0] - xq) / intv;
	double yq = y2.at[index_end - a - 1][0];
	return yq;
}

Matrix	gradient(Matrix _x, Matrix _y) {
	int mx = _x.rows;
	int my = _y.rows;
	int m = mx;
	double h = _x.at[1][0] - _x.at[0][0];
	Matrix df = createMat(_x.rows, 1);

	if (mx != my) {
		printf("Error: Dataset should have same point.\n");
	}
	if (m < 2) {
		printf("Error: Insufficient Dataset.\n");
	}
	else if (m == 2) {
		df.at[0][0] = (_y.at[1][0] - _y.at[0][0]) / h;
		df.at[1][0] = df.at[0][0];
	}
	else {
		df.at[0][0] = (4 * _y.at[1][0] - 3 * _y.at[0][0] - _y.at[2][0]) / (2 * h);
		for (int i = 1; i < m - 1; i++) {
			df.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / (2 * h);
		}
		df.at[m - 1][0] = (_y.at[m - 3][0] - 4 * _y.at[m - 2][0] + 3 * _y.at[m - 1][0]) / (2 * h);
	}
	return df;
}

void gradient1D(double x[], double y[], double dydx[], int m)
{
	double h = x[1] - x[0];
	if (m < 2) {
		printf("Error: Insufficient Dataset.\n");
	}
	else if (m == 2) {
		dydx[0] = (y[1] - y[0]) / h;
		dydx[1] = dydx[0];
	}
	else {
		dydx[0] = (4 * y[1] - 3 * y[0] -y[2]) / (2 * h);
		for (int i = 1; i < m - 1; i++) {
			dydx[i] = (y[i+1] - y[i-1]) / (2 * h);
		}
		dydx[m-1] = (y[m-3] - 4 * y[m - 2] + 3 * y[m - 1]) / (2 * h);
	}
}
Matrix	gradientFunc(double func(const double x), Matrix xin) {
	int m = xin.rows;
	Matrix y = createMat(m, 1);

	for (int i = 0; i < m; i++) {
		y.at[i][0] = func(xin.at[i][0]);
	}

	return gradient(xin, y);
}

void printArr(double x[], int m, const char* _name)
{
	printf("%s =\n", _name);
	for (int i = 0; i < m; i++)
	{
		printf("%15.6f\t", x[i]);
		printf("\n");
	}
	printf("\n");
}

double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), double _x0, double _tol)
{
	int k = 0;
	int Nmax = 100;
	double xk;
	double xk1;
	xk = _x0;
	double ep = 0;
	double hk;
	do {
		if (dfunc(xk) != 0) {                        // f'(xk) = 0 이 되면 Newton Raphson Method 사용 불가
			hk = -1 * func(xk) / dfunc(xk);    // hk 설정, xk = xk + hk
			xk1 = xk + hk;
			ep = fabs((xk1 - xk) / xk);
			k++;
			printf("Iteration:%d \t", k);
			printf("X(k): %f \t", xk);
			printf("Tolerance: %.10f\n", ep);
			xk = xk1;
		}
		else {
			printf("Wrong starting value");
			break;
		}

	} while (k<Nmax && ep>_tol);
	return xk;
}

double myFunc(const double x) {
	return  x * x * x;
}

double mydFunc(const double x) {
	return 3 * x * x;
}

Matrix	arr2Mat(double* _1Darray, int _rows, int _cols)
{
	Matrix Output = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
}