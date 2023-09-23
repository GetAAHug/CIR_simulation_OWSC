#include"initial_head.h"


double* Sphere_Bessel(double* x, int n, int flag)  //递归计算复数球bessel函数各阶数值。 flag = 1为第一类球bessel函数；flag = 2为第二类球bessel函数。
{


	complex<double> Piz (Pi,0);
	complex<double> X (x[0],x[1]);
	complex<double> c (n,0);
	complex<double> Bessel;
	//complex<double> Bessel_1;
	//complex<double> Bessel_2;
	double* result = new double[2];
	double* result_1 = new double[2];
	double* result_2 = new double[2];


	if (n == 0)
	{
		if (flag == 1)
		{
			Bessel = sqrt(complex<double>(2) / Piz / X)*sin(X);
			result[0] = real(Bessel);
			result[1] = imag(Bessel);
			return result;
		}
		else
		{
			Bessel = - sqrt(complex<double>(2) / Piz / X)*cos(X);
			result[0] = real(Bessel);
			result[1] = imag(Bessel);
			return result;
		}
	}
	else
	{
		if (n == -1)
		{
			if (flag == 1)
			{
				Bessel = sqrt(complex<double>(2) / Piz / X)*cos(X);
				result[0] = real(Bessel);
				result[1] = imag(Bessel);
				return result;
			}
			else
			{
				Bessel = sqrt(complex<double>(2) / Piz / X)*sin(X);
				result[0] = real(Bessel);
				result[1] = imag(Bessel);
				return result;
			}
		}
		else
		{
			result_1 = Sphere_Bessel(x, n - 1, flag);
			result_2 = Sphere_Bessel(x, n - 2, flag);
			complex<double> Bessel_1 (result_1[0], result_1[1]);
			complex<double> Bessel_2 (result_2[0], result_2[1]);
			Bessel = complex<double>(2.0 * n - 1) / X*Bessel_1 - Bessel_2;
			result[0] = real(Bessel);
			result[1] = imag(Bessel);
			return result;
		}
	}
		

}