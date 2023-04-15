#include "mod_matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void build_matrix_A(double dx, double dy, double** A, int Nx, int Ny){

  printf("build_A : begin");

  double b = -1/(dx*dx);
  double c = -1/(dy*dy);
  double a = 1 - 2*b - 2*c;

  // Initialisation de la matrice A

  for (int i = 0; i < Nx*Ny; i++){
    for (int j = 0; j < Nx*Ny; j++){
      A[i][j] = 0.0;
    }
  }

  // Remplissage de la matrice A

  for (int i = 1; i < Nx*Ny; i++){
    for (int j = 1; j < Nx*Ny; j++){
      if (i==j){
        A[i][j] = a;
        A[i][j-1] = b;
        A[i-1][j] = b;
        A[i][j-4] = c;
        A[i-4][j] = c;
      }
    }
  }

}
