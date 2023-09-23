//#include<cmath>
#include"initial_head.h"

double p_total(double cosalpha,double *atmos)
{ 
	double Ksray = atmos[0];
	double Ksmie = atmos[1];
	double Kamie = atmos[2];
	double Ks = Ksray + Ksmie;
	double g = atmos[3];
	double p_ray,p_mie,pfun;
	p_ray=(3*(1+3*gama+(1-gama)*pow(cosalpha,2.0)))/(16.0*Pi*(1+2*gama));
	p_mie=(1-pow(g,2.0))/(4.0*Pi)*(1/(pow((1+pow(g,2.0)-2*g*cosalpha),(3.0/2)))+f*0.5*(3*pow(cosalpha,2.0)-1)/pow((1+pow(g,2.0)),(3.0/2)));
	pfun=Ksray*p_ray/Ks+Ksmie*p_mie/Ks;
	return pfun;
}