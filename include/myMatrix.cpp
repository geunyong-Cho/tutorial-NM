/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [Cho Geunyong]
Created          : 26-03-2018
Modified         : 30-04-2021
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#include "myMatrix.h"
#include "myNM.h"


// Create Matrix with specified size
Matrix	createMat(int _rows, int _cols)
{
	// check matrix dimension
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out;
	// 1. Allocate row array first
	Out.at = (double**)malloc(sizeof(double*) * _rows);
	// 2. Then, allocate column 
	for (int i = 0; i < _rows; i++)
		Out.at[i] = (double*)malloc(sizeof(double) * _cols);
	// 3. Initialize row & column values of a matrix
	Out.rows = _rows;
	Out.cols = _cols;

	return Out;
}

// Create a matrix from a text file
Matrix	txt2Mat(std::string _filePath, std::string _fileName)
{
	std::ifstream file;
	std::string temp_string, objFile = _filePath + _fileName + ".txt";
	int temp_int = 0, nRows = 0;

	file.open(objFile);
	if (!file.is_open()) {
		printf("\n*********************************************");
		printf("\n  Could not access file: 'txt2Mat' function");
		printf("\n*********************************************\n");
		return createMat(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Out = createMat(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Out.at[i][j] = stof(temp_string);
		}
	file.close();

	return Out;
}

//basic

void  printMat(Matrix _A, const char* _name)
{
	printf("%s =\n", _name);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++)
			printf("%15.6f\t", _A.at[i][j]);
		printf("\n");
	}
	printf("\n");
}

Matrix copyMat(Matrix _A) {

	Matrix Out = createMat(_A.rows, _A.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[i][j] = _A.at[i][j];

	return Out;
}

Matrix transposeMat(Matrix _A) {

	Matrix Out = createMat(_A.cols, _A.rows);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[j][i] = _A.at[i][j];

	return Out;
}

Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}

Matrix multMat(Matrix _A, Matrix _B) {

	Matrix Out = createMat(_A.rows, _B.cols);

	if (_A.cols != _B.rows) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'multmat' function");
		printf("\n****************************************************\n");
	}
	else {

		for (int i = 0; i < _A.rows; i++)
		{
			for (int j = 0; j < _B.cols; j++)
			{
				Out.at[i][j] = 0;
				for (int k = 0; k < _A.cols; k++)
				{
					Out.at[i][j] = Out.at[i][j] + (_A.at[i][k] * _B.at[k][j]);
				}
			}
		}
		return Out;
	}
}

Matrix eye(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			if (i == j) {
				Out.at[i][j] = 1;
			}
			else {
				Out.at[i][j] = 0;
			}
	return Out;
}

Matrix eyevec(int _rows, int a) {
	Matrix Out = createMat(_rows, 1);
	for (int i = 0; i < _rows; i++)
	{
		if (i == a)
			Out.at[i][0] = 1;
		else
			Out.at[i][0] = 0;
	}
	return Out;
}

Matrix minus(Matrix _A)
{
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			_A.at[i][j] = -1 * _A.at[i][j];
	return _A;
}

void initMat(Matrix _A, double _val)
{
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			_A.at[i][j] = _val;
}

Matrix	zeros(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);

	initMat(Out, 0);

	return Out;
}

Matrix	ones(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);

	initMat(Out, 1);

	return Out;
}

// QR factorization
// norm of matrix
double normMat(Matrix _A)
{
	double norm = 0;
	for (int i = 0; i < _A.rows; i++) {
		norm += pow(_A.at[i][0], 2);
	}
	norm = sqrt(norm);
	return norm;
}

// 1 by 1 matrix divide
Matrix constdivideMat(Matrix _A, Matrix _B)
{
	Matrix Out = createMat(1, 1);
	Out.at[0][0] = _A.at[0][0] / _B.at[0][0];
	return Out;
}

// create vector c
Matrix Vector_c(Matrix _A, int a)
{
	Matrix Out = zeros(_A.rows, 1);
	for (int i = 0; i < _A.rows; i++) {
		Out.at[i][0] = _A.at[i][a];
	}
	return Out;
}

// create vector e
Matrix Vector_e(Matrix _A, int _k, double _ck)
{
	Matrix Out = createMat(_A.rows, 1);
	if (_k > _A.rows) {
		printf("Wrong input k");
	}
	for (int i = 0; i < _A.rows; i++)
		if (i != _k) {
			Out.at[i][0] = 0;
		}
		else {
			if (_ck >= 0) {
				Out.at[i][0] = 1;
			}
			else {
				Out.at[i][0] = -1;
			}
		}
	return Out;
}

// constant * matrix
Matrix multconstant(double _a, Matrix _A)
{
	Matrix Out = createMat(_A.rows, _A.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[i][j] = _a * _A.at[i][j];
	return Out;
}

Matrix Eigenvalue_withQR(Matrix _A, Matrix _Q, Matrix _R)
{
	Matrix Out = createMat(_A.rows, _A.cols);
	Matrix vecc = zeros(_A.rows, 1);
	Matrix vecv = copyMat(vecc);
	Matrix matH = zeros(_A.rows, _A.cols);
	Matrix matI = eye(_A.rows, _A.cols);
	Matrix matcoef = zeros(1, 1);
	Matrix coef = zeros(1, 1);
	initMat(matcoef, 2);
	Matrix vconst = zeros(1, 1);
	Matrix vece = copyMat(vecc);
	int a = 0;
	int Nmax = 1e+10;
	double tol = 1e-15;
	double ep;
	do {
		for (int k = 0; k < _R.rows - 1; k++) {
			vecc = Vector_c(_R, k);
			if (k >= 1) {
				for (int i = 0; i < k; i++) {
					vecc.at[i][0] = 0;
				}
			}
			vece = Vector_e(_A, k, vecc.at[k][0]);
			vecv = addMat(vecc, multconstant(normMat(vecc), vece));
			vconst = multMat(transposeMat(vecv), vecv);
			coef = constdivideMat(matcoef, vconst);
			matH = addMat(matI, minus(multconstant(coef.at[0][0], multMat(vecv, transposeMat(vecv)))));

			_Q = multMat(_Q, matH);
			_R = multMat(matH, _R);
		}
		ep = fabs((_A.at[0][0] - multMat(multMat(transposeMat(_Q), _A), _Q).at[0][0]) / _A.at[0][0]);
		_A = multMat(multMat(transposeMat(_Q), _A), _Q);
		_Q = eye(_A.rows, _A.cols);
		_R = copyMat(_A);
		a++;
	} while (a < Nmax && ep > tol);
	printMat(_A, "similar");
	printMat(eig_value(_A), "Eigen values");
	Out = eig_value(_A);
	return Out;
}

// find eigen values from similar matrix
Matrix eig_value(Matrix _A)
{
	Matrix Out = zeros(_A.rows, 1);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			if (i == j) {
				Out.at[_A.rows - i - 1][0] = _A.at[i][j];
			}
	return Out;
}

double cond(Matrix _A)
{
	double cond_A;
	Matrix trans = zeros(_A.cols, _A.rows);
	Matrix trans_eig = copyMat(trans);
	trans = multMat(transposeMat(_A), _A);
	Matrix Q = eye(trans.rows, trans.cols);	Matrix R = copyMat(trans);

	double a;
	double b;
	trans_eig = Eigenvalue_withQR(trans, Q, R);
	double sigma_max = fabs(sqrt(trans_eig.at[0][0]));
	double sigma_min = sigma_max;

	for (int i = 0; i < trans_eig.rows; i++) {
		a = fabs(sqrt(trans_eig.at[i][0]));
		if (sigma_max <= a) {
			sigma_max = a;
		}
	}
	for (int j = 0; j < trans_eig.rows; j++) {
		b = fabs(sqrt(trans_eig.at[j][0]));
		if (sigma_min == 0) {
			sigma_min = fabs(sqrt(trans_eig.at[j + 1][0]));
		}
		if (sigma_min >= b) {
			sigma_min = b;
		}
	}
	cond_A = sigma_max / sigma_min;
	return cond_A;
}

// swaping matrix values
void swapVal(Matrix _A, int a, int b) {
	double temp;

	for (int i = 0; i < _A.cols; i++) {
		temp = _A.at[a][i];
		_A.at[a][i] = _A.at[b][i];
		_A.at[b][i] = temp;
	}
}

// Apply back-substitution
Matrix	backSub(Matrix _A, Matrix _b)
{
	Matrix Out = createMat(_b.rows, 1);

	for (int i = _A.rows - 1; i >= 0; i--) {
		double temp = 0;
		for (int j = i + 1; j < _A.cols; j++)
			temp += _A.at[i][j] * Out.at[j][0];
		Out.at[i][0] = (_b.at[i][0] - temp) / _A.at[i][i];
	}
	return Out;
}

// Apply forward-substitution
Matrix	fwdSub(Matrix _A, Matrix _b)
{
	Matrix Out = createMat(_b.rows, 1);

	for (int i = 0; i < _A.rows; i++) {
		double temp = 0;
		for (int j = 0; j < i; j++)
			temp += _A.at[i][j] * Out.at[j][0];
		Out.at[i][0] = (_b.at[i][0] - temp) / _A.at[i][i];
	}
	return Out;
}

Matrix gaussElimU(Matrix _A, Matrix _b, Matrix _U, Matrix _d)
{
	int M = _A.rows;
	int N = _A.cols;
	Matrix m = createMat(M, 1);
	_U = _A;
	_d = _b;
	if (M == N)
	{
		for (int k = 0; k < M - 1; k++)
		{
			for (int i = k + 1; i < M; i++)
			{
				if (_U.at[k][k] != 0)
				{
					m.at[i][0] = _U.at[i][k] / _U.at[k][k];
				}
				else {
					printf("Zero-division error!!");
					break;
				}

				for (int j = k; j < N; j++)
				{
					_U.at[i][j] = _U.at[i][j] - m.at[i][0] * _U.at[k][j];
				}
				_d.at[i][0] = _d.at[i][0] - m.at[i][0] * _d.at[k][0];
			}
		}
	}
	else {
		printf("Rows and Columns of Matrix A should be same");
	}

	return _U;
}

Matrix gaussElimd(Matrix _A, Matrix _b, Matrix _U, Matrix _d)
{
	int M = _A.rows;
	int N = _A.cols;
	Matrix m = createMat(M, 1);
	_U = _A;
	_d = _b;
	if (M == N)
	{
		for (int k = 0; k < M - 1; k++)
		{
			for (int i = k + 1; i < M; i++)
			{
				if (_U.at[k][k] != 0)
				{
					m.at[i][0] = _U.at[i][k] / _U.at[k][k];
				}
				else {
					printf("Zero-division error!!");
					break;
				}

				for (int j = k; j < N; j++)
				{
					_U.at[i][j] = _U.at[i][j] - m.at[i][0] * _U.at[k][j];
				}
				_d.at[i][0] = _d.at[i][0] - m.at[i][0] * _d.at[k][0];
			}
		}
	}
	else {
		printf("Rows and Columns of Matrix A should be same");
	}

	return _d;
}

void gaussjordan(Matrix _A, Matrix _b, Matrix _U, Matrix _d)
{
	int m = _A.rows;
	int n = _A.cols;
	double temp1;
	double temp2;
	for (int k = 0; k < m; k++) {
		temp1 = _U.at[k][k];
		_d.at[k][0] = _d.at[k][0] / temp1;
		for (int j = k; j < n; j++) {
			_U.at[k][j] = _U.at[k][j] / temp1;
		}
		for (int i = 0; i < m; i++) {
			temp2 = _U.at[i][k];
			if (i != k) {
				_d.at[i][0] = _d.at[i][0] - temp2 * _d.at[k][0];
				for (int j = k; j < n; j++) {
					_U.at[i][j] = _U.at[i][j] - temp2 * _U.at[k][j];
				}
			}
		}
	}
	printMat(_U, "Matrix U");
	printMat(_d, "Vector d");
}

void LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P)
{
	int M = _A.rows;
	int N = _A.cols;
	double max = 0;
	int piv_row;
	double a;
	double temp = 0;
	double sp;
	Matrix vec = createMat(_A.rows, 1);
	Matrix m = createMat(M, 1);
	Matrix I = eye(M, N);
	if (M == N)
	{
		for (int k = 0; k < M; k++)
		{
			// pivoting process
			//piv_row = find_pivot(_U, k);
			for (int i = k; i < M; i++) {
				a = abs(_U.at[i][k]);
				for (int j = k; j < N; j++) {
					if (max < abs(_U.at[i][j])) {
						max = abs(_U.at[i][j]);
					}
					sp = a / max;
					if (temp < sp) {
						piv_row = i;
					}
					temp = sp;
				}
			}
			if (piv_row != k) {
				swapVal(_L, k, piv_row); swapVal(_U, k, piv_row); swapVal(_P, k, piv_row);
				/*		printf("At k = %d\n\n", k);
						printMat(_P, "P");
						printMat(_L, "L");
						printMat(_U, "U");*/
			}
			for (int i = k + 1; i < M; i++) {
				if (_U.at[k][k] != 0)
				{
					m.at[i][0] = _U.at[i][k] / _U.at[k][k];
					// matrix L
					_L.at[i][k] = m.at[i][0];
				}
				else {
					printf("Zero-division error!!");
					break;
				}
				for (int j = k; j < N; j++)
				{
					_U.at[i][j] = _U.at[i][j] - m.at[i][0] * _U.at[k][j];
				}
				// matrix U
				_U.at[i][k] = 0;
			}
		}
	}
	else {
		printf("Rows and Columns of Matrix A should be same");
	}
}

// solving linear equation
Matrix solveLU(Matrix _L, Matrix _U, Matrix _b, Matrix _P)
{
	Matrix a = createMat(_L.rows, 1);
	a = multMat(_P, _b);
	Matrix _y = createMat(_L.rows, 1);
	Matrix _x = createMat(_L.rows, 1);
	_y = fwdSub(_L, a);
	_x = backSub(_U, _y);
	return _x;
}

// inverse Matrix
Matrix inv(Matrix _A, Matrix _Ainv)
{
	Matrix L = zeros(_A.rows, _A.cols);
	Matrix U = copyMat(_A);
	Matrix P = eye(_A.rows, _A.cols);
	Matrix b = createMat(_A.rows, 1);
	Matrix _x = zeros(_A.rows, 1);
	Matrix I = eye(_A.rows, _A.cols);

	LUdecomp(_A, L, U, P);
	L = addMat(L, I);

	for (int i = 0; i < _A.cols; i++) {
		b = eyevec(_A.rows, i);
		_x = solveLU(L, U, b, P);
		for (int j = 0; j < _A.cols; j++) {
			_Ainv.at[j][i] = _x.at[j][0];
		}
	}
	return _Ainv;
}

