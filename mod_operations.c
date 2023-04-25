#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mod_operations.h"

double* vector_sum(double* u, double* v, int n)
// Sum of u+v = x
{
    double* x = (double*)malloc(n*sizeof(double));

    for (int i = 0; i < n; i++) {
        x[i] = u[i] + v[i];
    }

    return x;
}


double* vector_substract(double* u, double* v, int n)
// Substraction of u-v = x
{
    double* x = (double*)malloc(n*sizeof(double));

    for (int i = 0; i < n; i++) {
        x[i] = u[i] - v[i];
    }

    return x;
}


double* vector_coef(double* u, double c, int n)
// Apply coefficient to c*u = x
{
    double* x = (double*)malloc(n*sizeof(double));

    for (int i = 0; i < n; i++) {
        x[i] = c * u[i];
    }

    return x;
}


double vector_scalar(double* u, double* v, int n)
// Dot product of u and v
{
    double x = 0;

    for (int i = 0; i < n; i++) {
        x += u[i] * v[i];
    }

    return x;
}


double* matvect_product(double* x, int Nx, int Ny, double DeltaT, double DeltaX, double DeltaY)
// Matrix vector product of A (matrix) and x (vector)
{
    double* b = (double*)malloc(Nx*Ny*sizeof(double));

    double C, Cx, Cy;

    C = 1 + 2*DeltaT/(DeltaX*DeltaX) + 2*DeltaT/(DeltaY*DeltaY);
    Cx = -DeltaT/(DeltaX*DeltaX);
    Cy = -DeltaT/(DeltaY*DeltaY);

    for (int i = 0; i < Nx*Ny; i++) {
        b[i] = C * x[i];
        if (i%Nx != Nx-1) {
            b[i] += Cx * x[i+1];
        }
        if (i%Nx != 0) {
            b[i] += Cx * x[i-1];
        }
        if (i/Nx != Ny-1) {
            b[i] += Cy * x[i+Nx];
        }
        if (i/Nx != 0) {
            b[i] += Cy * x[i-Nx];
        }
    }

    return b;
}
