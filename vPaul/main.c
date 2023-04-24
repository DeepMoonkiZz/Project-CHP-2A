#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mod_function.h"
#include "mod_matrix.h"
#include "mod_gradient.h"
#include "mod_display.h"

int main() {

  // Lecture des paramètre dans le fichier parameter.dat
  // ---------------------------------------------------------------------------

  FILE *file_parameter;
  int Nx, Ny, kmax;
  double Lx, Ly, D, dt, Tmax, eps;
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
  sscanf(val_file, "%d\n%d\n%lf\n%lf\n%lf\n%lf\n%d\n%lf\n%lf\n", &Nx, &Ny, &Lx, &Ly, &D, &dt, &kmax, &Tmax, &eps);
  fclose(file_parameter);

  // Initialisation
  // ---------------------------------------------------------------------------

  // Calcul des pas de maillage

  double dx = Lx/(double)(Nx + 1);
  double dy = Ly/(double)(Ny + 1);

  // Allocation et initialisation de la mémoire des vecteurs principaux

  double * u_0 = (double*)malloc((Nx*Ny+8)*sizeof(double));
  double * u_1 = (double*)malloc((Nx*Ny+8)*sizeof(double));

  for (int i = 0; i < Nx*Ny+8; i++){
    u_0[i] = 0.0;
    u_1[i] = 0.0;
  }

  // Début de la boucle en utilisant le conjugate_gradiant
  // ---------------------------------------------------------------------------

  double t = 0;

  while(t < Tmax){
    conjugate_gradiant(u_1,u_0,kmax,D,t,dt,dx,dy,eps,Nx,Ny);
    for (int i = 0; i < Nx*Ny+8; i++){
      u_0[i] = u_1[i];
      u_1[i] = 0.0;
    }
    t += dt;
    printf("%lf\n", t);
  }
  display_u(u_0,Nx,Ny);

  // Libération de la mémoire allouée

  free(u_0);
  free(u_1);

  return 0;

}
