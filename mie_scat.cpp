#include"initial_head.h"

double* mie_scat(double lambda, double radius, double density, bool flag)//计算大气mie氏散射相关参数(Ks_mie,Ka_mie,g)。flag = 0乡村，flag = 1城市。
{
	double* M = new double[2];

	M = refractive_index(lambda, flag);
	double* mie_p = new double[3];

	double* mie_parameters = new double[3];

	mie_parameters = mie_cal(M, lambda, radius);
	mie_p[0] = mie_parameters[0] * Pi * density * pow(radius, 2) / (1e6);
	mie_p[1] = (mie_parameters[1] - mie_parameters[0]) * Pi*density * pow(radius, 2) / (1e6);
	mie_p[2] = mie_parameters[2] / mie_parameters[0];
	delete[]mie_parameters;

	return mie_p;
}