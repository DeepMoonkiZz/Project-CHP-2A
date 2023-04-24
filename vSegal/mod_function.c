#include <math.h>

#include "mod_function.h"

#define PI acos(-1)


double func(double x, double y,double Lx, double Ly, double t, int f)
{
  if (f==1) {
    return 2*(y - y*y + x - x*x);
  }
  else if (f==2) {
    return sin(x) + cos(y);
  }
  else {
    return exp(-(x-Lx/2)*(x-Lx/2))*exp(-(y-Ly/2)*(y-Ly/2))*cos(t*PI/2);
  }
}


double g(double x, double y, int f)
{
  if (f==1) {
    return 0;
  }
  else if (f==2) {
    return sin(x) + cos(y);
  }
  else {
    return 0;
  }
}


double h(double x, double y, int f)
{  
  if (f==1) {
    return 0;
  }
  else if (f==2) {
    return sin(x) + cos(y);
  }
  else {
    return 1;
  }
}