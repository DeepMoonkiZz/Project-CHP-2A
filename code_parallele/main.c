#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include <openmpi/mpi.h>
#include "/home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/fonction/charge.h"

#include "mod_gradient.h"
#include "mod_operations.h"
#include "mod_scheme.h"
#include "mod_function.h"
#include "/home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/parameter_and_display/mod_display.h"
#include "/home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/parameter_and_display/mod_parameter.h"


int main(int argc, char *argv[])
{
    // ---------------------------------------------------------------------------

    // Lecture et définitions des paramètre dans le fichier parameter.dat

    // ---------------------------------------------------------------------------

    // Enregistrer l'heure de début
    clock_t start = clock();

    // Initialisation de la parallelisation

    int rank, nproc;    
    int iBeg, iEnd;
    iBeg = 0;
    iEnd = 0;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Création de la structure

    struct data data;

    Read_parameter(&data);

    charge(rank, data.Nx*data.Ny, nproc, &iBeg, &iEnd);

    data.rank = rank, data.npart = iEnd-iBeg+1, data.nproc = nproc, data.iBeg = iBeg, data.iEnd = iEnd;

    // Définitions des tableau de stockage de u ainsi que du second membre et du temps
    
    double  *u_exact_part, *u_part, *b_part, *u, *u_exact;

    u_part = (double*)malloc(data.npart*sizeof(double));
    b_part = (double*)malloc(data.npart*sizeof(double));
    u_exact_part = (double*)malloc(data.npart*sizeof(double));

    u = (double*)malloc(data.Nx*data.Ny*sizeof(double));
    u_exact = (double*)malloc(data.Nx*data.Ny*sizeof(double));

    FILE *file_display = 0;
    // Définition du fichier de sortie de l'erreur
    if (data.rank == 0) {
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
            sprintf(filename, "Solutions/Extension_2/Erreur/error.dat");
        }
        else {
            printf("Erreur dans le numero de fonction choisit.\n");
            exit(1);
        }
        file_display = fopen(filename, "w");
    }

    // ---------------------------------------------------------------------------

    // Execution du programme de résolution du problème

    // ---------------------------------------------------------------------------

    // Initialisation du probleme
    for (int i = data.iBeg; i <= data.iEnd; i++) {
        u_part[i-iBeg] = sol_init((i%data.Nx)*data.DeltaX, (i/data.Nx)*data.DeltaY, data);
        b_part[i-iBeg] = 0.;
    }

    // Boucle en temps pour les itérations

    int p = 0;
    double t = 0, error = 0;
    
    while (data.Tmax > t + pow(10,-10)) {
        t += data.DeltaT;
        
        // Compute u
        Build_vect_b(b_part, u_part, t, data);
        gradient_conjugate(u_part, b_part, data);

        // Compute u exact
        Build_u_exact(u_exact_part, t, data);

        // Send u_part and u_exact_part to proc 0
        Build_comp_vect(u, u_part, data);
        Build_comp_vect(u_exact, u_exact_part, data);

        error = Calculate_error(u_exact, u, data);
        if (data.rank == 0) {
            // Compute error and save in file
            fprintf(file_display, "%lf\t%lf\n", t, error);

            // Save u and u_exact in file for code_sequentieleach dt
            if (data.function == 3) {
                display_u(u, data, p, 1);
                display_u(u_exact, data, p, 2);
                p++;
            }
        }
    }
    


    // Save u and u_exact in file
    if (data.rank == 0) {
        fclose(file_display);
        if (data.function != 3) {
            display_u(u, data, p, 1);
            display_u(u_exact, data, p, 2);
        }
    }

    // Deallocate everything
    free(u_part), free(b_part), free(u_exact), free(u);
    MPI_Finalize();

    // Enregistrer l'heure de fin
    clock_t end = clock();

    // Calculer le temps écoulé en secondes (en utilisant le rapport CLOCKS_PER_SEC)
    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;

    // Afficher le temps écoulé
    printf("Le code a mis %.2f secondes à s'exécuter.\n", elapsed_time);

    return 0;
}