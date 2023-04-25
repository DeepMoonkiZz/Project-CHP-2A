#include <stdio.h>
#include <stdlib.h>

#include "mod_display.h"


void display_u(double * u, int Nx, int Ny, double DeltaX, double DeltaY, int p){

  // Le but de cette fonction est de construire le fichier qui sera lu
  // par Gnuplot.

  FILE *file_display = 0;
  char filename[100];

  // ------------- PARTIE A MODULER -------------

  sprintf(filename, "Solutions/Instationnaire/sol_%d.dat", p);
  file_display = fopen(filename, "w");

  // --------------------------------------------

  int k = 0;

  for(int j = 0; j < Ny; j++){
    for(int i = 0; i < Nx; i++){
      fprintf(file_display, "%lf\t%lf\t%lf\n", i*DeltaX, j*DeltaY, u[k]);
      k++;
    }
  }
  fclose(file_display);

}
