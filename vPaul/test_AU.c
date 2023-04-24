#include <stdlib.h>
#include <stdio.h>
#include "mod_function.h"

double f_1(double x, double y){

  return 2*(y - y*y + x - x*x);
}

int main(){

  double dx  = 0.1;
  double dy = 0.1;
  double dt = 0.01;
  double D = 1.;
  double b = 1./(dx*dx);
  double c = 1./(dy*dy);
  double a = 1. - 2*b - 2*c;
  int Nx = 3;
  int Ny = 2;
  int Imax = Nx*Ny;
  int l;

  double * v = malloc(Imax*sizeof(double));
  double * u = malloc((Imax+8)*sizeof(double));
  double * u_1 = malloc(Imax*sizeof(double));
  double * F = malloc(Imax*sizeof(double));
  double * vect_b = malloc(Imax*sizeof(double));
  double * Aub = malloc(Imax*sizeof(double));

  // NB : pour éviter les problèmes de bords, on va décaler les indices et
  // remplir le vecteur b avec une "marge de 0.0" à droite et à gauche.


  for (int i = 0; i < 5; i++){
    u[i] = 0.0;
  }

  for (int i = Imax ; i < Imax+5; i++){
    u[i] = 0.0;
  }

  for (int i = 0; i < Imax; i++){
    u[i+4] = 2.;
    u_1[i] = 0.0;
    v[i] = 0.0;
  }

  for (int k = 0; k < Imax; k++){

      if((k+1)%(Nx) != 0){
          v[k] = c*u[k] + b*u[k+3] + a*u[k+4] + b*u[k+5] + c*u[k+8];
      }
      else{
          v[k] = c*u[k] + b*u[k+3] + a*u[k+4] + c*u[k+8];
          k++;
          v[k] = c*u[k] + a*u[k+4] + b*u[k+5] + c*u[k+8];
      }
  }

  printf("Vecteur u = \n");
  for (int i = 0; i < Imax+8; i++){
    printf("u[%d] = %lf\n", i, u[i]);
  }

  printf("-----\n");
  printf("\n");

  printf("Vecteur Au = \n");
  for (int i = 0; i < Imax; i++){
    printf("Au[%d] = %lf\n", i, v[i]);
  }

  for (int j = 1; j < Ny+1; j++){
    for (int i = 1; i < Nx+1; i++){
      l = Nx*(j-1) + (i-1);
      F[l] = f_1(i*dx,j*dy);
    }
  }

  printf("-----\n");
  printf("\n");

  printf("Vecteur F = \n");
  for (int i = 0; i < Imax; i++){
    printf("F[%d] = %lf\n", i, F[i]);
  }

  for (int i = 0; i < Imax; i ++){
    vect_b[i] = (1./dt)*u[i+4] + (1./D)*F[i];
  }

  printf("-----\n");
  printf("\n");

  printf("Vecteur vect_b = \n");
  for (int i = 0; i < Imax; i++){
    printf("vect_b[%d] = %lf\n", i, vect_b[i]);
  }

  for (int i = 0; i < Imax; i ++){
    Aub[i] = vect_b[i] - v[i];
  }

  printf("-----\n");
  printf("\n");
  printf("Vecteur Au - b = \n");
  for (int i = 0; i < Imax; i++){
    printf("Aub[%d] = %lf\n", i, Aub[i]);
  }

  free(u);
  free(u_1);
  free(v);
  free(vect_b);
  free(F);

  return 0;
}
