#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <openmpi/mpi.h>
#include "/home/segal/Documents/MatMeca/S8/CHP/fonction/charge.h"

#include "mod_gradient.h"
#include "mod_operations.h"
#include "mod_scheme.h"
#include "mod_function.h"
#include "/home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/parameter_and_display/mod_display.h"
#include "/home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/parameter_and_display/mod_parameter.h"


int main(int argc, char ** argv)
{
    // ---------------------------------------------------------------------------

    // Lecture et définitions des paramètre dans le fichier parameter.dat

    // ---------------------------------------------------------------------------

    // Création de la structure

    struct data data; // = {.rank = rank, .nproc = nproc, .iBeg = iBeg, .iEnd = iEnd};

    Read_parameter(&data);

    // Définitions des tableau de stockage de u ainsi que du second membre et du temps

    double *u, *u_exact, *b;
    u = (double*)malloc(data.Nx*data.Ny*sizeof(double));
    u_exact = (double*)malloc(data.Nx*data.Ny*sizeof(double));
    b = (double*)malloc(data.Nx*data.Ny*sizeof(double));

    double t = 0, error = 0;

    // Définition du fichier de sortie de l'erreur

    FILE *file_display = 0;
    char filename[100];
    if (data.function==1) {
        sprintf(filename, "Solutions/Stationnaire_1/Erreur/error.dat");
    }
    else if (data.function==2) {
        sprintf(filename, "Solutions/Stationnaire_2/Erreur/error.dat");
    }
    else if (data.function==3) {
        sprintf(filename, "Solutions/Instationnaire/Erreur/error.dat");   
    }
    else if (data.function==4) {
        sprintf(filename, "Solutions/Extension_1/Erreur/error.dat");
    }
    else if (data.function==5) {
        sprintf(filename, "Solutions/Extension_1/Erreur/error.dat");
    }
    else {
        printf("Erreur dans le numero de fonction choisit.\n");
        exit(1);
    }
    file_display = fopen(filename, "w");

    // ---------------------------------------------------------------------------

    // Execution du programme de résolution du problème

    // ---------------------------------------------------------------------------

    // Initialisation du probleme
    for (int j = 0; j < data.Ny; j++) {
        for (int i = 0; i < data.Nx; i++) {
            u[data.Nx*j + i] = sol_init(i*data.DeltaX, j*data.DeltaY, data);
            b[data.Nx*j + i] = 0.;
        }
    }

    // Boucle en temps pour les itérations

    int p = 0;

    while (t<data.Tmax) {
        t += data.DeltaT;

        // Compute u
        Build_vect_b(b, u, t, data);
        gradient_conjugate(u, b, data);

        // Compute u exact
        Build_u_exact(u_exact, t, data);

        // Compute error and save in file
        error = Calculate_error(u_exact, u, data.Nx*data.Ny);
        fprintf(file_display, "%lf\t%lf\n", t, error);

        // Save u and u_exact in file for code_sequentieleach dt
        if (data.function == 3) {
            display_u(u, data, p, 1);
            display_u(u_exact, data, p, 2);
            p++;
        }
    }

    // Save u and u_exact in file
    if (data.function != 3) {
        display_u(u, data, p, 1);
        display_u(u_exact, data, p, 2);
    }

    // Deallocate everything
    free(u_exact), free(u), free(b);
    fclose(file_display);

    return 0;
}