/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University
Author           : [ Cho Geunyong ]
Created          : 30-04-2021
Modified         : 15-05-2021
Language/ver     : C++ in MSVS2019

Description      : NM_Assignment6_21600676_Á¶±Ù¿ë.cpp
-------------------------------------------------------------------------------*/
#define Assignment	6    	// enter your assignment number
#define eval		0		// set 0

#include "../../include/myNM.h"
#include <math.h>

int main(int argc, char* argv[])
{
	/*	 [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif
	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*==========================================================================*/

	double t_arr[] = { 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
	double pos_arr[] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59 };
	double vel_arr[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double acc_arr[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double xin_arr[] = { 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0 };
	double x0 = 3;
	double tol = 1e-5;
	int M = 21;
	Matrix t = arr2Mat(t_arr, M, 1);
	Matrix pos = arr2Mat(pos_arr, M, 1);

	/*==========================================================================*/
	/*					      Numerical method algorithm     	                */
	/*==========================================================================*/
	printf("\n**************************************************");
	printf("\n|                     PART 1.                    |");
	printf("\n**************************************************\n");

	Matrix vel = gradient(t, pos);
	Matrix acc = gradient(t, vel);
	printMat(t, "t");
	printMat(pos, "pos");
	printMat(vel, "vel_mat");
	printMat(acc, "acc_mat");

	gradient1D(t_arr, pos_arr, vel_arr, M);
	printArr(vel_arr, M, "vel_arr");
	gradient1D(t_arr, vel_arr, acc_arr, M);
	printArr(acc_arr, M, "acc_arr");

	// PART 2
	printf("\n**************************************************");
	printf("\n|                     PART 2.                    |");
	printf("\n**************************************************\n");

	Matrix xin = arr2Mat(xin_arr, M, 1);
	Matrix dydx = gradientFunc(myFunc, xin);
	printMat(xin, "xin");
	printMat(dydx, "dydx");

	double xk = newtonRaphsonFunc(myFunc, mydFunc, x0, tol);

	system("pause");
	return 0;
}


