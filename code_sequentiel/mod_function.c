#include <math.h>

#include "mod_function.h"

#define PI acos(-1)


double func_exact(double x, double y, double t, struct data dt)
{
  if (dt.function==1) {
    return x*(1-x)*y*(1-y);
  }
  else if (dt.function==2) {
    return sin(x) + cos(y);
  }
  else if (dt.function==3) {
    return exp(-(x-dt.Lx/2)*(x-dt.Lx/2))*exp(-(y-dt.Ly/2)*(y-dt.Ly/2))*cos(t*PI/2);
  }
  else {
    return cos(3.*PI*x)*(y-(dt.ymax+dt.ymin)/4.)/(t+1.);
  }
}


double func(double x, double y, double t, struct data dt)
{
  if (dt.function==1) {
    return 2*(y - y*y + x - x*x);
  }
  else if (dt.function==2) {
    return sin(x) + cos(y);
  }
  else if (dt.function==3) {
    return exp(-(x-dt.Lx/2)*(x-dt.Lx/2))*exp(-(y-dt.Ly/2)*(y-dt.Ly/2))*cos(t*PI/2);
  }
  else {
    return -cos(3*PI*x)/pow(t+1.,2)*(y-(dt.ymin+dt.ymax)/4)+9*dt.D*pow(PI,2)*cos(3*PI*x)*(y-(dt.ymax+dt.ymin)/4)/(t+1.);
  }
}


double g(double x, double y, struct data dt)
{
  if (dt.function==1) {
    return 0;
  }
  else if (dt.function==2) {
    return sin(x) + cos(y);
  }
  else {
    return 0;
  }
}


double h(double x, double y, struct data dt)
{
  if (dt.function==1) {
    return 0;
  }
  else if (dt.function==2) {
    return sin(x) + cos(y);
  }
  else {
    return 1;
  }
}


double h_left(double y, double t, struct data dt)
{
  return func_exact(dt.xmin-dt.DeltaX, y, t, dt);
  // return -3*PI*sin(3.*PI*dt.xmin)*(y-(dt.ymin+dt.ymax)/4.)/(t+1.0);
}

double g_right(double y, double t, struct data dt)
{
  return func_exact(dt.xmax+dt.DeltaX, y, t, dt);
}

double h_top(double x, double t, struct data dt)
{
  return func_exact(x, dt.ymax+dt.DeltaY, t, dt);
  // return cos(3.*PI*x)/(t+1.0);
}

double g_bottom(double x, double t, struct data dt)
{
  return func_exact(x, dt.ymin-dt.DeltaY, t, dt);
}

double sol_init(double x, double y, struct data dt)
{
  if (dt.function==1) {
    return 0.;
  }
  else if (dt.function==2) {
    return 0.;
  }
  else if (dt.function==3) {
    return 0.;
  }
  else {
    return func_exact(x, y, 0, dt);
  }
}