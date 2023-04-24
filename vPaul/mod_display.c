#include "mod_display.h"
#include <stdio.h>

void display_u(double * u, int Nx, int Ny){

  // Le but de cette fonction est de construire le fichier qui sera lu
  // par Gnuplot.

  FILE *file_display = fopen("solution.dat","w");
  int k = 4;

  for(int j = 0; j < Ny; j++){
    for(int i = 0; i < Nx; i++){
      fprintf(file_display, "%lf\t", u[k]);
      k++;
    }
    fprintf(file_display, "\n");
  }

}
