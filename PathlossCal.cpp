#include"initial_head.h"

double PathlossCal(double* atmos, double* geo, double* parameters, double Sr, int M, int n, int N)
{
	double sita_T = geo[0];
	double sita_R = geo[1];
	double beta_T = geo[2];
	double beta_R = geo[3];
	double fai_T = geo[4];
	double fai_R = geo[5];
	double r = geo[6];

	double Ks = atmos[0] + atmos[1];
	double Ke = atmos[0] + atmos[1] + atmos[2];

	double epsilon = parameters[0];
	double delta = parameters[1];

	double PD_i[5] = { 0 };
	double u_T[3] = { cos(fai_T)*sin(sita_T),sin(fai_T)*sin(sita_T),cos(sita_T) };
	double u_R[3] = { cos(fai_R)*sin(sita_R),sin(fai_R)*sin(sita_R),cos(sita_R) };
	double r_R = sqrt(Sr / Pi);
	const double c = 2.997046e8;

	if ((-u_T[1] >= cos(beta_T / 2)) && (u_R[1] >= cos(beta_R / 2)))
	{
		PD_i[0] = exp(-Ke*r)*(1 - r*pow((pow(r, 2) + pow(r_R, 2)), (-1.0 / 2))) / (1 - cos(beta_T / 2));
	}
	int i;
	double location[3], u[3], u_new[3], location_new[3], location_minus_new[3];
	double d_total, p, alpha, psi, temp, d, Pi_R, PD_cur;
	bool flag;

	for (int k = 1; k <= N; k++)
	{
		location[0] = 0;
		location[1] = r;
		location[2] = 0;
		u[0] = u_T[0];
		u[1] = u_T[1];
		u[2] = u_T[2];
		d_total = 0;
		p = 1 - PD_i[0];
		i = 0;
		while ((i < n) && (p >= epsilon))
		{
			i = i + 1;
			flag = false;
			while (!flag)
			{
				if (i == 1)
				{
					alpha = acos(1 - rand() / (RAND_MAX + 0.0)*(1 - cos(beta_T / 2)));
					psi = 2 * Pi*rand() / (RAND_MAX + 0.0);
				}
				else
				{
					alpha = scat_angle(M, delta, atmos);
					psi = 2 * Pi*rand() / (RAND_MAX + 0.0);
				}
				temp = sqrt(1 - pow(u[2], 2));
				if (!(temp == 0.0))
				{
					u_new[0] = sin(alpha)*(u[0] * u[2] * cos(psi) - u[1] * sin(psi)) / temp + u[0] * cos(alpha);
					u_new[1] = sin(alpha)*(u[1] * u[2] * cos(psi) + u[0] * sin(psi)) / temp + u[1] * cos(alpha);
					u_new[2] = -sin(alpha)*cos(psi)*temp + u[2] * cos(alpha);
				}
				else
				{
					u_new[0] = sin(alpha)*cos(psi);
					u_new[1] = sin(alpha)*sin(psi);
					u_new[2] = u[2] * cos(alpha);
				}
				d = -log(1 - rand() / (RAND_MAX + 0.0)) / Ke;
				location_new[0] = location[0] + d*u_new[0];
				location_new[1] = location[1] + d*u_new[1];
				location_new[2] = location[2] + d*u_new[2];
				location_minus_new[0] = location[0] - location_new[0];
				location_minus_new[1] = location[1] - location_new[1];
				location_minus_new[2] = location[2] - location_new[2];
				double* det = new double[3];
				if ((((norm(cross(location, location_new, det))) / norm(location_minus_new)) > r_R) || (dot(location, u_new) >= 0) || (d < norm(location)))
					flag = true;
				delete[]det;
			}
			location[0] = location_new[0];
			location[1] = location_new[1];
			location[2] = location_new[2];
			u[0] = u_new[0];
			u[1] = u_new[1];
			u[2] = u_new[2];
			d_total = d_total + d;
			Pi_R = min(1.0, 2 * Pi*(1 - norm(location)*pow((pow(norm(location), 2) + pow(r_R, 2)), -1.0 / 2))
				*p_total(-dot(location, u) / norm(location), atmos))*exp(-Ke*norm(location));
			if ((dot(location, u_R) / norm(location)) >= cos(beta_R / 2))
			{
				PD_cur = p * Ks / Ke*Pi_R;
				PD_i[i] = PD_i[i] + PD_cur;
			}
			p = p * Ks / Ke*(1 - Pi_R);
		}
	}


	double PD_overall = PD_i[0], PL_overall;
	for (int m = 1; m < n + 1; m++)
	{
		PD_i[m] = PD_i[m] / N;
		PD_overall = PD_overall + PD_i[m];
	}

	PL_overall = 10 * log10(1.0 / PD_overall);
	return PL_overall;
}
