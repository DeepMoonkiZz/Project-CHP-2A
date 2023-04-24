#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mod_gradient.h"
#include "mod_operations.h"
#include "mod_scheme.h"
#include "mod_display.h"


int main(int argc, char ** argv) 
{
    // ---------------------------------------------------------------------------

    // Lecture et définitions des paramètre dans le fichier parameter.dat 
    
    // ---------------------------------------------------------------------------

    FILE *file_parameter;
    int Nx, Ny, kmax, f;
    double Lx, Ly, D, DeltaT, Tmax, eps;
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
    sscanf(val_file, "%d\n%d\n%lf\n%lf\n%lf\n%lf\n%d\n%lf\n%lf\n%d\n", &Nx, &Ny, &Lx, &Ly, &D, &DeltaT, &kmax, &Tmax, &eps, &f);
    fclose(file_parameter);

    // Définitions de la variable de temps ainsi que les pas d'espace

    double t = 0, DeltaX = Lx/(Nx-1), DeltaY = Ly/(Ny-1);

    // Définitions des tableau de stockage de u ainsi que du second membre

    double *u, *b; 
    u = (double*)malloc(Nx*Ny*sizeof(double));
    b = (double*)malloc(Nx*Ny*sizeof(double));


    // ---------------------------------------------------------------------------

    // Execution du programme de résolution du problème

    // ---------------------------------------------------------------------------

    // Initialisation du probleme 
    for (int i=0; i<Nx*Ny; i++) {
        u[i] = 0.;
    }

    // Boucle en temps pour les itérations
    while (t<Tmax) {
        t += DeltaT;
        b = Build_vect_b(u, t, Nx, Ny, DeltaT, DeltaX, DeltaY, Lx, Ly, D, f);
        u = gradient_conjugate(u, b, eps, kmax, Nx, Ny, DeltaT, DeltaX, DeltaY);
        free(b);
    }

    display_u(u, Nx, Ny, DeltaX, DeltaY);


    return 0;
}