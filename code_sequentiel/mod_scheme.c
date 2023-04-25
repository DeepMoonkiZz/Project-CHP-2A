#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mod_function.h"
#include "mod_scheme.h"


void Build_vect_b(double* b, double* u, double t, int Nx, int Ny, double DeltaT, double DeltaX, double DeltaY, double Lx, double Ly, double D, int f)
{
    double Cx = DeltaT*D / (DeltaX*DeltaX);
    double Cy = DeltaT*D / (DeltaY*DeltaY);

    memcpy(b, u, Nx*Ny*sizeof(double));

    // Définition des boundary conditions

    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            b[Nx*j + i] += DeltaT * func(i*DeltaX, j*DeltaY, Lx, Ly, t, f);
        }
    }

    for (int i = 0; i < Nx; i++) {
        b[i] += Cx * g(i*DeltaX, 0, f);
        b[Nx*(Ny-1) + i] += Cx * g(i*DeltaX, Ly, f);
    }

    for (int j = 0; j < Ny; j++) {
        b[Nx*j] += Cy * h(0, j*DeltaY, f);
        b[Nx*j + Ny-1] += Cy * h(Lx, j*DeltaX, f);
    }
}