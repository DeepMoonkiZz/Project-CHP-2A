#include <stdio.h>
#include <stdlib.h>

#include "mod_display.h"


void display_u(double * u, struct data dt, int p, int type_sol){

  // Le but de cette fonction est de construire le fichier qui sera lu
  // par Gnuplot.

  FILE *file_display = 0;
  char filename[100];

  // ------------- OUVERTURE DU FICHIER EN FONCTION DE LA FONCTION -------------

  if (type_sol == 1) {
    if (dt.function == 1) {
      sprintf(filename, "Solutions/Stationnaire_1/Numerique/sol_%d.dat", p);
    }

    else if (dt.function == 2) {
      sprintf(filename, "Solutions/Stationnaire_2/Numerique/sol_%d.dat", p);
    }

    else if (dt.function == 3) {
      sprintf(filename, "Solutions/Instationnaire/Numerique/sol_%d.dat", p);
    }
    else if (dt.function == 4) {
      sprintf(filename, "Solutions/Extension_1/Numerique/sol_%d.dat", p);
    }
    else if (dt.function == 5) {
      sprintf(filename, "Solutions/Extension_2/Numerique/sol_%d.dat", p);
    }
  }

  else if (type_sol == 2) {
    if (dt.function == 1) {
      sprintf(filename, "Solutions/Stationnaire_1/Exact/sol_%d.dat", p);
    }

    else if (dt.function == 2) {
      sprintf(filename, "Solutions/Stationnaire_2/Exact/sol_%d.dat", p);
    }

    else if (dt.function == 3) {
      sprintf(filename, "Solutions/Instationnaire/Exact/sol_%d.dat", p);
    }    
    else if (dt.function == 4) {
      sprintf(filename, "Solutions/Extension_1/Exact/sol_%d.dat", p);
    }    
    else if (dt.function == 5) {
      sprintf(filename, "Solutions/Extension_2/Exact/sol_%d.dat", p);
    }    
  }

  file_display = fopen(filename, "w");

  // ---------------------------------------------------------------------------

  int k = 0;

  for(int j = 0; j < dt.Ny; j++){
    for(int i = 0; i < dt.Nx; i++){
      fprintf(file_display, "%lf\t%lf\t%lf\n", dt.xmin + i*dt.DeltaX, dt.ymin + j*dt.DeltaY, u[k]);
      k++;
    }
  }
  fclose(file_display);

}
