#include "mod_function.h"
#include <stdio.h>
#include <math.h>
#define PI 3.1415926536

double f_1(double x, double y){

  return 2*(y - y*y + x - x*x);
}

double f_2(double x, double y){

  return sin(x) + cos(y);
}

double f_3(double x, double y,double Lx, double Ly, double t){

  return exp(-(x-Lx/2)*(x-Lx/2))*exp(-(y-Ly/2)*(y-Ly/2))*cos(t*PI/2);
}

double g_null(){

  return 0;
}

double g(double x, double y){

  return sin(x) + cos(y);
}

double h_const(int val){

  if(val==0){
    return 0;
  }
  else if(val==1){
    return 1;
  }
  else{
    printf("ERROR : h_const(val) mal configurée");
    return 0;
  }

}

double h(double x, double y){

  return sin(x) + cos(y);
}
