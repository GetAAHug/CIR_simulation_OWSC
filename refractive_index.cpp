#include"initial_head.h"

double* refractive_index(double lambda, bool flag)//计算大气的复折射率。flag = 0是乡村，flag = 1是城市。
{
	double Ws_real, Ws_imag;
	double Dl_real, Dl_imag;
	double Sl_real, Sl_imag;
	double* M = new double[2];

	//Calculate water soluble
	if (lambda < 0.7)
	{
		Ws_real = 1.53;
	}
	else
	{
		if (lambda < 1.2)
		{
			Ws_real = 1.52;
		}
		else
		{
			Ws_real = 1.51;
		}
		
	}

	if (lambda < 0.3371)
	{
		Ws_imag = 3.8168 * pow(lambda, 3) + 0.7374 * pow(lambda, 2) - 1.7139 * lambda + 0.3527;
	}
	else
	{
		if (lambda < 0.52)
		{
			Ws_imag = 0.005;
		}
		else
		{
			if (lambda < 0.632)
			{
				Ws_imag = 0.006;
			}
			else
			{
				Ws_imag = -0.1684 * pow(lambda, 5) + 0.9966 * pow(lambda, 4) - 2.2945 * pow(lambda, 3) + 2.5492 * pow(lambda, 2) - 1.3372 * lambda + 0.2701;
			}
		}
	}
	//Calculate dust like
	if (lambda < 0.7)
	{
		Dl_real = 1.53;
	}
	else
	{
		if (lambda < 1.06)
		{
			Dl_real = 1.52;
		}
		else
		{
			Dl_real = -0.0089 * pow(lambda, 2) - 0.229 * lambda + 1.7727;
		}
	}

	if (lambda < 0.3)
	{
		Dl_imag = 3.6 * pow(lambda, 2) - 2.42 * lambda + 0.41;
	}
	else
	{
		Dl_imag = 0.008;
	}
	//Calculate soot like
	if (lambda < 0.3371)
	{
		Sl_real = -160.0396 * pow(lambda, 3) + 118.4297 * pow(lambda, 2) - 26.4073 * lambda + 3.3246;
	}
	else
	{
		if (lambda < 1.06)
		{
			Sl_real = 1.75;
		}
		else
		{
			Sl_real = 1.75 + 0.042 * (lambda - 1.06);
		}
	}

	if (lambda < 0.25)
	{
		Sl_imag = 0.35 + 2 * (lambda - 0.2);
	}
	else
	{
		Sl_imag = 0.45;
	}

	complex<double> Ws(Ws_real, Ws_imag);
	complex<double> Dl(Dl_real, Dl_imag);
	complex<double> Sl(Sl_real, Sl_imag);
	complex<double> m;

	//Calculate the refractive index
	if (flag == false)
	{
		m = complex<double>(0.3) * Ws + complex<double>(0.7) * Dl;
	}
	else
	{
		if (flag == true)
		{
			m = complex<double>(0.8) * (complex<double>(0.3) * Ws + complex<double>(0.7) * Dl) + complex<double>(0.2) * Sl;
		}
	}
		


	M[0] = m.real();
	M[1] = m.imag();

	return M;
}




