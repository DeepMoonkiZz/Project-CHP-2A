#include "mod_parameter.h"

void Read_parameter(struct data *dt)
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
    }

    //Récupération des valeures à l'aide de la chaîne de charactères val_file

    fgets(val_file, 100, file_parameter);
    fgets(val_file, 100, file_parameter);
    sscanf(val_file, "%d\n%d\n%lf\n%lf\n%lf\n%lf\n%d\n%lf\n%lf\n%d\n", &Nx, &Ny, &Lx, &Ly, &D, &DeltaT, &kmax, &Tmax, &eps, &f);
    fclose(file_parameter);

    dt->function = f;
    dt->Kmax = kmax;
    dt->eps = eps;
    dt->Tmax = Tmax;

    dt->Nx = Nx;
    dt->Ny = Ny;
    dt->Lx = Lx;
    dt->Ly = Ly;

    dt->D = D;
    
    dt->DeltaT = DeltaT; 
    dt->DeltaX = Lx/(Nx-1);
    dt->DeltaY = Ly/(Ny-1);
}