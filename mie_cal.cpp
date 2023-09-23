#include"initial_head.h"

double* mie_cal(double* M, double lambda, double radius)   //返回(Qsca,Qext,Qcos)的地址
{
	double x = 2 * Pi * radius / lambda;
	const int Nmax = 15;      //球bessel级数的前 Nmax 项求和来近似
	double* Mx = new double[2];
	double* X = new double[2];
	X[0] = x;
	X[1] = 0;
	complex<double> m(M[0], M[1]);

	Mx[0] = M[0] * x;
	Mx[1] = M[1] * x;
	complex<double> mx(Mx[0], Mx[1]);

	/*存放各类bessel函数值*/
	double(* bm)[Nmax] = new double[2][Nmax];
	double* b1 = new double[Nmax];
	double* b2 = new double[Nmax];
	double(* h)[Nmax] = new double[2][Nmax];

	/*存放bessel函数求导中间变量*/
	double(* bm_s)[Nmax] = new double[2][Nmax];
	double* b1_s = new double[Nmax];
	double* b2_s = new double[Nmax];
	double(* h_s)[Nmax] = new double[2][Nmax];

	/*存放各类bessel函数的导数值*/
	double(* bm_d)[Nmax] = new double[2][Nmax];
	double* b_d = new double[Nmax];
	double(* h_d)[Nmax] = new double[2][Nmax];

	for (int i = 0; i < Nmax; i++)
	{
		double* temp1 = new double[2];
		double* temp2 = new double[2];
		double* temp3 = new double[2];

		temp1 = Sphere_Bessel(Mx, i + 1, 1);
		temp2 = Sphere_Bessel(X, i + 1, 1);
		temp3 = Sphere_Bessel(X, i + 1, 2);

		complex<double> bmc(temp1[0], temp1[1]);
		complex<double> sqz = sqrt(complex<double>(0.5*Pi) / mx);
		complex<double> bmz = bmc*sqz;
		bm[0][i] = bmz.real();
		bm[1][i] = bmz.imag();

		b1[i] = temp2[0] * sqrt(0.5*Pi / x);

		b2[i] = temp3[0] * sqrt(0.5*Pi / x);

		h[0][i] = b1[i];
		h[1][i] = b2[i];

		delete[]temp1;
		delete[]temp2;
		delete[]temp3;
	}

	for (int i = 0; i < Nmax; i++)
	{
		if (i == 0)
		{
			complex<double> bz0 = sin(mx) / mx;
			bm_s[0][i] = bz0.real();
			bm_s[1][i] = bz0.imag();
			b1_s[i] = sin(x) / x;
			b2_s[i] = -cos(x) / x;
			h_s[0][i] = b1_s[0];
			h_s[1][i] = b2_s[0];
		}
		else
		{
			bm_s[0][i] = bm[0][i - 1];
			bm_s[1][i] = bm[1][i - 1];
			b1_s[i] = b1[i - 1];
			b2_s[i] = b2[i - 1];
			h_s[0][i] = b1_s[i];
			h_s[1][i] = b2_s[i];
		}
	}
	
	for (int i = 0 ; i < Nmax; i++)
	{
		b_d[i] = x*b1_s[i] - (i + 1)*b1[i];

		complex<double> bmz(bm[0][i], bm[1][i]);
		complex<double> bm_sz(bm_s[0][i], bm_s[1][i]);
		complex<double> bm_dz = mx*bm_sz - complex<double>(i + 1)*bmz;
		bm_d[0][i] = bm_dz.real();
		bm_d[1][i] = bm_dz.imag();

		complex<double> hz(h[0][i], h[1][i]);
		complex<double> h_sz(h_s[0][i], h_s[1][i]);
		complex<double> h_dz = complex<double>(x)*h_sz - complex<double>(i + 1)*hz;
		h_d[0][i] = h_dz.real();
		h_d[1][i] = h_dz.imag();
	}
	delete[]bm_s;
	delete[]b1_s;
	delete[]b2_s;
	delete[]h_s;


	double(* parameters)[2][Nmax] = new double[2][2][Nmax];

	for (int i = 0; i < Nmax; i++)
	{
		complex<double> bmz(bm[0][i], bm[1][i]);
		complex<double> bm_dz(bm_d[0][i], bm_d[1][i]);
		complex<double> bz(b1[i], 0);
		complex<double> b_dz(b_d[i], 0);
		complex<double> hz(h[0][i], h[1][i]);
		complex<double> h_dz(h_d[0][i], h_d[1][i]);

		complex<double> a = (pow(m, 2)*bmz*b_dz - bz*bm_dz) / (pow(m, 2)*bmz*h_dz - hz*bm_dz);
		complex<double> b = (bmz*b_dz - bz*bm_dz) / (bmz*h_dz - hz*bm_dz);
		//complex<double> c = (bz*h_dz - hz*b_dz) / (bmz*h_dz - hz*bm_dz);
		//complex<double> d = m*(bz*h_dz - hz*b_dz) / (pow(m, 2)*bmz*h_dz - hz*bm_dz);

		parameters[0][0][i] = a.real();
		parameters[0][1][i] = a.imag();
		parameters[1][0][i] = b.real();
		parameters[1][1][i] = b.imag();
		//parameters[2][0][i] = c.real();
		//parameters[2][1][i] = c.imag();
		//parameters[3][0][i] = d.real();
		//parameters[3][1][i] = d.imag();
	}

	delete[]bm;
	delete[]b1;
	delete[]b2;
	delete[]h;
	delete[]bm_d;
	delete[]b_d;
	delete[]h_d;

	double* result = new double[3];
	for (int i = 0; i < 3; i++)
	{
		result[i] = 0;
	}

	for (int j = 0; j < 2; j++)
	{
		for (int i = 0; i < Nmax; i++)
		{
			result[0] = result[0] + (pow(parameters[j][0][i], 2) + pow(parameters[j][1][i], 2))*(2 * i + 3) * 2 / pow(x, 2);
		}
	}

	for (int j = 0; j < 2; j++)
	{
		for (int i = 0; i < Nmax; i++)
		{
			result[1] = result[1] + parameters[j][0][i] * (2 * i + 3) * 2 / pow(x, 2);
		}
	}

	for (int i = 0; i < Nmax - 1; i++)
	{
		complex<double> a_now(parameters[0][0][i], parameters[0][1][i]);
		complex<double> a_next(parameters[0][0][i + 1], parameters[0][1][i + 1]);
		complex<double> b_now(parameters[1][0][i], parameters[1][1][i]);
		complex<double> b_next(parameters[1][0][i + 1], parameters[1][1][i + 1]);
		complex<double> sum = a_now*conj(a_next) + b_now*conj(b_next);
		complex<double> prod = a_now*conj(b_now);

		result[2] = result[2] + (sum.real()*(i + 1)*(i + 3) / (i + 2) + prod.real()*(2 * i + 3) / (i + 1) / (i + 2)) * 4 / pow(x, 2);
	}
	delete[]parameters;

	return result;

}