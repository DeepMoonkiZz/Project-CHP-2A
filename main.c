#include <stdio.h>
#include <stdlib.h>

#include "mod_function.h"
#include "mod_matrix.h"

int main() {

  // Lecture des paramètre dans le fichier parameter.dat
  // ---------------------------------------------------------------------------

  FILE *file_parameter;
  int Nx, Ny;
  double Lx, Ly, D, dt;
  char val_file[100];

  file_parameter = fopen("parameter.dat","r");

  // Verification de l'existence du fichier

  if(file_parameter == NULL) {
    printf("Erreur : Impossible d'ouvrir le fichier parameter.dat.\n");
    return 1;
  }

  //Récupération des valeures à l'aide de la chaîne de charactères val_file

  fgets(val_file, 100, file_parameter);
  fgets(val_file, 100, file_parameter);
  sscanf(val_file, "%d\n%d\n%lf\n%lf\n%lf\n%lf\n", &Nx, &Ny, &Lx, &Ly, &D, &dt);

  fclose(file_parameter);

  // Initialisation
  // ---------------------------------------------------------------------------

  // Calcul des pas de maillage

  double dx = Lx/(double)(Nx + 1);
  double dy = Ly/(double)(Ny + 1);

  // Allocation et remplissage de la matrice A

  double **A = malloc((Nx*Ny)*sizeof(double*));
  for (int i = 0; i < Nx; i++) {
    A[i] = malloc((Nx*Ny)*sizeof(double));
  }

  return 0;

}
