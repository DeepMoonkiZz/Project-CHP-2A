#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mod_operations.h"
#include "mod_gradient.h"
#include "mod_function.h"


double* Build_vect_b(double* u, double t, int Nx, int Ny, double DeltaT, double DeltaX, double DeltaY, double Lx, double Ly, double D)
{
    double* b = (double*)malloc(Nx*Ny*sizeof(double));
    double Cx = DeltaT*D / (DeltaX*DeltaX);
    double Cy = DeltaT*D / (DeltaY*DeltaY);

    b = u;

    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            b[Nx*j + i] += DeltaT * f_3(i*DeltaX, j*DeltaY, Lx, Ly, t);
        }
    }

    for (int i=0; i<Nx; i++) {
        b[i] += Cx * g(i*DeltaX, 0);
        b[Nx*(Ny-1) + i] += Cx * g(i*DeltaX, Ly);
    }

    for (int j=0; j<Ny; j++) {
        b[Nx*j] += Cy * h(0, j*DeltaY);
        b[Nx*j + Ny-1] += Cy * h(Lx, j*DeltaX);
    }

    return b;
}