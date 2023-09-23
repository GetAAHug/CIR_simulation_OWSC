#include"initial_head.h"

double ray_scat(double lambda)
{
	double Ns = 2.54743e25;
	double m = (5791817 / (238.0185 - pow(lambda, -2)) + 167909 / (57.362 - pow(lambda, -2)))*1e-8 + 1;
	double dep = 0.036;
	double Ks_ray = 24 * pow(Pi, 3) / pow(lambda, 4) / Ns * 1e24 * pow(pow(m, 2) - 1, 2) / pow(pow(m, 2) + 2, 2)*(6 + 3 * dep) / (6 - 7 * dep);
	return Ks_ray;
}