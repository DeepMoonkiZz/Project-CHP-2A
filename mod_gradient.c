#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mod_operations.h"
#include "mod_gradient.h"

double* gradient_conjugate(double* x, double* b, double eps, double kmax, int Nx, int Ny, double DeltaT, double DeltaX, double DeltaY)
{
    int k = 0, n = Nx*Ny;

    double alpha, gamma, beta;
    double* r = (double*)malloc(n*sizeof(double));
    double* p = (double*)malloc(n*sizeof(double));
    double* z = (double*)malloc(n*sizeof(double));
    double* xplus = (double*)malloc(n*sizeof(double));
    double* rplus = (double*)malloc(n*sizeof(double));

    r = vector_substract(b, matvect_product(x, Nx, Ny, DeltaT, DeltaX, DeltaY), n);
    p = r;
    beta = sqrt(vector_scalar(r, r, n));

    while (beta > eps && k < kmax) {

        z = matvect_product(p, Nx, Ny, DeltaT, DeltaX, DeltaY); // z = Ap
        alpha = vector_scalar(r, r, n) / vector_scalar(z, p, n); // alpha

        xplus = vector_sum(x, vector_coef(p, alpha, n), n); // x += alpha*x
        rplus = vector_substract(r, vector_coef(z, alpha, n), n); // r -= alpha*r

        gamma = vector_scalar(rplus, rplus, n) / vector_scalar(r, r, n); // gamma
        p = vector_sum(rplus, vector_coef(p, gamma, n), n); // p = r + gamma*p

        beta = sqrt(vector_scalar(rplus, rplus, n)); // beta = norme(r)

        r = rplus;
        x = xplus;
        k++;

    }
    if (k > kmax) {
        printf("Tolérence non atteinte norme=%f", beta);
    }

    return x;
}
