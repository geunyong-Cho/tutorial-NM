/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [Cho Geunyong]
Created          : 26-03-2018
Modified         : 30-04-2021
Language/ver     : C++ in MSVS2019

Description      : myMatrix.h
----------------------------------------------------------------*/

#ifndef		_MY_MATRIX_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>

typedef struct {
	double** at;
	int rows, cols;
}Matrix;

// Create Matrix with specified size
extern	Matrix createMat(int _rows, int _cols);

// Create a matrix from a text file
// basic
extern	Matrix txt2Mat(std::string _filePath, std::string _fileName);

extern	void printMat(Matrix _A, const char* _name);

extern Matrix copyMat(Matrix _A);

extern Matrix transposeMat(Matrix _A);

extern Matrix addMat(Matrix _A, Matrix _B);

extern Matrix multMat(Matrix _A, Matrix _B);

extern Matrix eye(int _rows, int _cols);

extern Matrix minus(Matrix _A);

extern void initMat(Matrix _A, double _val);

extern Matrix zeros(int _rows, int _cols);

extern Matrix ones(int _rows, int _cols);

// QR factorization
extern double normMat(Matrix _A);

extern Matrix constdivideMat(Matrix _A, Matrix _B);

extern Matrix Vector_c(Matrix _A, int a);

extern Matrix Vector_e(Matrix _A, int _k, double _ck);

extern Matrix multconstant(double _a, Matrix _A);

extern Matrix Eigenvalue_withQR(Matrix _A, Matrix _Q, Matrix _R);

extern Matrix eig_value(Matrix _A);

extern double cond(Matrix _A);

// LU decomposition
extern void swapVal(Matrix _A, int a, int b);

extern	Matrix	backSub(Matrix _A, Matrix _b);

extern	Matrix  fwdSub(Matrix _A, Matrix _b);

extern Matrix gaussElimU(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

extern Matrix gaussElimd(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

void gaussjordan(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

extern void LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P);

extern Matrix solveLU(Matrix _L, Matrix _U, Matrix _b, Matrix _P);

extern Matrix inv(Matrix _A, Matrix _Ainv);

#endif