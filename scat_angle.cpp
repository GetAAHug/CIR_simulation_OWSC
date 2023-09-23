#include"initial_head.h"

double scat_angle(int M, double delta, double* atmos)
{
	double Ksray = atmos[0];
	double Ksmie = atmos[1];
	double Kamie = atmos[2];
	double Ks = Ksray + Ksmie;
	double g = atmos[3];

	double* C = new double[6];
	C[0] = -Ksray / Ks * 3.0 / 8 * (1.0 + 3 * gama) / (1.0 + 2 * gama) + Ksmie / Ks * f / 4.0 * (1 - pow(g, 2.0)) / pow((1 + pow(g, 2.0)), (3.0 / 2));
	C[1] = -Ksray / Ks / 8.0 * (1.0 - gama) / (1.0 + 2 * gama) - Ksmie / Ks / 4.0 * f * (1 - pow(g, 2.0)) / pow((1 + pow(g, 2.0)), (3.0 / 2));
	C[2] = -Ksmie / Ks * (1.0 - pow(g, 2.0)) / (2.0 * g);
	C[3] = 1 + pow(g, 2.0);
	C[4] = 2 * g;
	C[5] = 1.0 - Ksray / Ks * (1.0 / 2) + Ksmie / Ks * (1.0 - g) / (2 * g);

	double alpha;
	double u = rand() / (RAND_MAX + 0.0);
	double v = 0;
	bool flag = false;
	int k = 0;
	double tmp = Fv(v, C) - u;
	while ((k < M) && (!flag))
	{
		k = k + 1;
		v = v - tmp / (2.0 * Pi * p_total(v, atmos));
		v = min(1.0, max(-1.0, v));
		tmp = Fv(v, C) - u;
		if (fabs(tmp) <= delta)
			flag = true;
	}
	if (!flag)
	{
		v = 0;
		double step = 1;
		while (step > delta)
		{
			step = step / 2.0;
			if (Fv(v, C) <= u)
				v = v + step;
			else
				v = v - step;
		}
	}
	alpha = acos(v);
	delete[]C;
	return alpha;
}