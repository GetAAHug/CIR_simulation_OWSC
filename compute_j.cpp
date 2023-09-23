#include"initial_head.h"
int compute_j(double t,double T_i[],int length)
{
	int j = 0;
	for (int i = 0; i < length - 1; i++)
	{
		if ((t >= T_i[i]) && (t < T_i[i + 1]))
		{
			j = i + 1;
		}
	}
	return j;
}