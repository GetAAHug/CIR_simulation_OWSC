//#include<cmath>
#include"initial_head.h"

double Fv(double v, double* C)
{
	double F;
	F = 1 - (C[0] * v + C[1] * pow(v, 3.0) + C[2] * pow((C[3] - C[4] * v), (-1.0 / 2)) + C[5]);
	return F;
}

