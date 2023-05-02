#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <openmpi/mpi.h>
#include "/home/segal/Documents/MatMeca/S8/CHP/fonction/charge.h"

#include "mod_gradient.h"
#include "mod_operations.h"
#include "mod_scheme.h"
#include "mod_display.h"


int main(int argc, char ** argv)
{
    // ---------------------------------------------------------------------------

    // Lecture et définitions des paramètre dans le fichier parameter.dat

    // ---------------------------------------------------------------------------

    // Intialisation de MPI et de la fonction charge

/*    int rank, nproc;    
    int iBeg, iEnd;

    iBeg = 0;
    iEnd = 0;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
*/

    // Création de la structure

    struct data data; // = {.rank = rank, .nproc = nproc, .iBeg = iBeg, .iEnd = iEnd};

    Read_parameter(&data);

    // Définitions des tableau de stockage de u ainsi que du second membre et du temps

    double *u, *b;
    u = (double*)malloc(data.Nx*data.Ny*sizeof(double));
    b = (double*)malloc(data.Nx*data.Ny*sizeof(double));

    double t;

    // ---------------------------------------------------------------------------

    // Execution du programme de résolution du problème

    // ---------------------------------------------------------------------------

    // Initialisation du probleme
    for (int i = 0; i < data.Nx*data.Ny; i++) {
        u[i] = 0.;
        b[i] = 0.;
    }

    // Boucle en temps pour les itérations

    int p = 0;

    while (t<data.Tmax) {
        t += data.DeltaT;
        Build_vect_b(b, u, t, data);
        gradient_conjugate(u, b, data);
        // display_u(u, Nx, Ny, DeltaX, DeltaY, p);
        // p++;
    }

    display_u(u, data.Nx, data.Ny, data.DeltaX, data.DeltaY, p);


    return 0;
}