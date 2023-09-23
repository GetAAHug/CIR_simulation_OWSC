#include<iostream>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<fstream>
#include<algorithm>
#include<complex>
#include<math.h>
#include<omp.h>
#include<windows.h>
#include<stdlib.h>

using namespace std;

const double Pi = 3.141592653589793;
const double f = 0.5;
const double gama = 0.017;

double Fv(double,double*);
double scat_angle(int, double, double*);
double p_total(double, double*);
double norm(double*);
double dot(double*, double*);
double* cross(double*, double*, double*);
double* mie_cal(double*, double, double);
double* mie_scat(double, double, double, bool);
double* refractive_index(double, bool);
double* Sphere_Bessel(double*, int, int);
double ray_scat(double);
double PathlossCal(double*, double*, double*, double, int, int, int);

int compute_j(double, double[], int);
