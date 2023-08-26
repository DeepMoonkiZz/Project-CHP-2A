#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mod_operations.h"


void vector_sum(double* x, double* u, double* v, int n)
// Sum of u+v = x
{
    for (int i = 0; i < n; i++) {
        x[i] = u[i] + v[i];
    }
}


void vector_substract(double* x, double* u, double* v, int n)
// Substraction of u-v = x
{
    for (int i = 0; i < n; i++) {
        x[i] = u[i] - v[i];
    }
}


void vector_coef(double* x, double* u, double c, int n)
// Apply coefficient to c*u = x
{
    for (int i = 0; i < n; i++) {
        x[i] = c * u[i];
    }
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


void matvect_product(double* b, double* x, int Nx, int Ny, double DeltaT, double DeltaX, double DeltaY, double D, int f)
// Matrix vector product of A (laplacian matrix) and x (vector) giving the result b
{
    double C, Cx, Cy;

    Cx = -DeltaT*D/(DeltaX*DeltaX);
    Cy = -DeltaT*D/(DeltaY*DeltaY);
    C = 1 - 2*Cx - 2*Cy;

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
/*
    if (f==4 || f==5) {
        for (int i = 0; i < Nx; i++) {
            b[Nx*(Ny-1) + i] -= Cy * x[Nx*(Ny-1) + i];
        }
        for (int j = 0; j < Ny; j++) {
            b[Nx*j] -= Cx * x[Nx*j];
        }
    }
*/
}


double Calculate_error(double* u_exact, double* u, int n)
{
    double sum=0;
    for (int i=0; i<n; i++) {
        sum += (u_exact[i]-u[i])*(u_exact[i]-u[i]);
    }
    sum = sqrt(sum);
    sum = sum/(1.0*n);
    
    return sum;
}