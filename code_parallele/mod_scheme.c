#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mod_function.h"
#include "mod_scheme.h"


void Build_vect_b(double* b, double* u, double t, struct data dt)
{
    double Cx = dt.DeltaT*dt.D / (dt.DeltaX*dt.DeltaX);
    double Cy = dt.DeltaT*dt.D / (dt.DeltaY*dt.DeltaY);

    // Remplissage du vecteur b avec les composantes du vecteur U

    memcpy(b, u, dt.Nx*dt.Ny*sizeof(double));

    // Remplissage du vecteur b avec les composantes du vecteur F

    for (int j = 0; j < dt.Ny; j++) {
        for (int i = 0; i < dt.Nx; i++) {
            b[dt.Nx*j + i] += dt.DeltaT * func(i*dt.DeltaX, j*dt.DeltaY, dt.Lx, dt.Ly, t, dt.function);
        }
    }

    // Remplissage du vecteur b avec les composantes des conditions limites

    for (int i = 0; i < dt.Nx; i++) {
        b[i] += Cx * g(i*dt.DeltaX, 0, dt.function);
        b[dt.Nx*(dt.Ny-1) + i] += Cx * g(i*dt.DeltaX, dt.Ly, dt.function);
    }

    for (int j = 0; j < dt.Ny; j++) {
        b[dt.Nx*j] += Cy * h(0, j*dt.DeltaY, dt.function);
        b[dt.Nx*j + dt.Ny-1] += Cy * h(dt.Lx, j*dt.DeltaX, dt.function);
    }
}