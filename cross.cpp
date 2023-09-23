#include"initial_head.h"
double* cross(double*a,double*b,double* c)
{
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[2]*b[0]-a[0]*b[2];
	c[2]=a[0]*b[1]-a[1]*b[0];
	return c;
}