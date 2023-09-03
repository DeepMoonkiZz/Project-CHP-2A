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
            b[dt.Nx*j + i] += dt.DeltaT * func(dt.xmin + i*dt.DeltaX, dt.ymin + j*dt.DeltaY, t, dt);
        }
    }

    // Remplissage du vecteur b avec les composantes des conditions limites

    if (dt.function==4 || dt.function==5) {
        // Up and down side
        for (int i = 0; i < dt.Nx; i++) {
            b[i] += Cy * g_bottom(dt.xmin + i*dt.DeltaX, t, dt);
            b[dt.Nx*(dt.Ny-1) + i] -= dt.DeltaY * Cy * h_top(dt.xmin + i*dt.DeltaX, t, dt);
        }
        // Left and right side
        for (int j = 0; j < dt.Ny; j++) {
            b[dt.Nx*j] -= dt.DeltaX * Cx * h_left(dt.ymin + j*dt.DeltaY, t, dt); 
            b[dt.Nx*j + dt.Nx-1] += Cx * g_right(dt.ymin + j*dt.DeltaY, t, dt);
        }
    }
    else {
        // Up and down side
        for (int i = 0; i < dt.Nx; i++) {
            b[i] += Cy * g(i*dt.DeltaX, dt.ymin - dt.DeltaY, dt);
            b[dt.Nx*(dt.Ny-1) + i] += Cy * g(i*dt.DeltaX, dt.ymax + dt.DeltaY, dt);
        }
        // Left and right side
        for (int j = 0; j < dt.Ny; j++) {
            b[dt.Nx*j] += Cx * h(dt.xmin - dt.DeltaX, j*dt.DeltaY, dt);
            b[dt.Nx*j + dt.Nx-1] += Cx * h(dt.xmax + dt.DeltaX, j*dt.DeltaX, dt);
        }
    }
}


void Build_u_exact(double* u_exact, double t, struct data dt)
{
    for (int j = 0; j < dt.Ny; j++) {
        for (int i = 0; i < dt.Nx; i++) {
            u_exact[dt.Nx*j + i] = func_exact(dt.xmin+i*dt.DeltaX, dt.ymin+j*dt.DeltaY, t, dt);
        }
    }
}