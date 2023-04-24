#include "mod_display.h"
#include <stdio.h>

void display_u(double * u, int Nx, int Ny, double DeltaX, double DeltaY){

  // Le but de cette fonction est de construire le fichier qui sera lu
  // par Gnuplot.

  FILE *file_display = fopen("solution.dat","w");
  int k = 0;

  for(int j = 0; j < Ny; j++){
    for(int i = 0; i < Nx; i++){
      fprintf(file_display, "%lf\t%lf\t%lf\n", i*DeltaX, j*DeltaY, u[k]);
      k++;
    }
  }

}
