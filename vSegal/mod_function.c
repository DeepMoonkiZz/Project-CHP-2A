#include <math.h>

#include "mod_function.h"

#define PI acos(-1)


double f_1(double x, double y)
{
  return 2*(y - y*y + x - x*x);
}


double f_2(double x, double y)
{
  return sin(x) + cos(y);
}


double f_3(double x, double y,double Lx, double Ly, double t)
{
  return exp(-(x-Lx/2)*(x-Lx/2))*exp(-(y-Ly/2)*(y-Ly/2))*cos(t*PI/2);
}


double g(double x, double y)
{
  return 0;
  // return sin(x) + cos(y);
}


double h(double x, double y)
{
  // return 0;
  // return sin(x) + cos(y);
  return 1;
}
