/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [Cho Geunyong]
Created          : 26-03-2018
Modified         : 30-04-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H
#define		PI		3.14159265358979323846264338327950288419716939937510582

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myMatrix.h"

// New function add
void gradient2(double xin[]);

// BisectionNL
extern double bisection(double _a0, double _b0, double _tol);

// Newton Raphson
extern double newtonRaphson(double _x0, double _tol);

// Secant Method
extern double secant(double _x1, double _x2, double _tol);

//Original Function
extern double original_func(double _x);

// 1 diff Function
extern double ass1_func(double _x);

// 2 diff Function
extern double ass1_dfunc(double _x);

// Hybrid Method
extern double hybrid_method(double _a0, double _b0, double _x0, double _tol);

// Bonus Newton
extern double bonus_newtonRaphson(double _x0, double _tol);

// Bonus Function
extern double bonus_func(double _x);

// Bonus 1diff Function
extern double bonus_dfunc(double _x);

extern double func1(double _x, double _y);

extern double func2(double _x, double _y);

extern double dfunc1x(double _x, double _y);

extern double dfunc1y(double _x, double _y);

extern double dfunc2x(double _x, double _y);

extern double dfunc2y(double _x, double _y);

extern Matrix FX(double _x, double _y);

extern Matrix Jac(double _x, double _y);

extern void System_non_linear(double _x, double _y, double _tol);

extern Matrix linearFit(Matrix _x, Matrix _y);

extern double linearinterp(Matrix _x, Matrix _y, int xq);

extern Matrix gradient(Matrix _x, Matrix _y);

extern Matrix gradientFunc(double func(const double x), Matrix xin);

extern void gradient1D(double x[], double y[], double dydx[], int m);

extern void printArr(double x[], int m, const char* _name);

extern double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), double _x0, double _tol);

extern double myFunc(const double x);

extern double mydFunc(const double x);

Matrix arr2Mat(double* _1Darray, int _rows, int _cols);
#endif

