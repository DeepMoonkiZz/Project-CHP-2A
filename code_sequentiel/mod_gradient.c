#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mod_operations.h"
#include "mod_gradient.h"

void gradient_conjugate(double* u, double* b, double eps, double kmax, int Nx, int Ny, double DeltaT, double DeltaX, double DeltaY)
{
    int k = 0, n = Nx*Ny;

    double alpha, gamma, beta;
    double* r = (double*)malloc(n*sizeof(double));
    double* p = (double*)malloc(n*sizeof(double));
    double* z = (double*)malloc(n*sizeof(double));
    double* uplus = (double*)malloc(n*sizeof(double));
    double* rplus = (double*)malloc(n*sizeof(double));
    double* alphap = (double*)malloc(n*sizeof(double));
    double* alphaz = (double*)malloc(n*sizeof(double));
    double* gammap = (double*)malloc(n*sizeof(double));
    
    matvect_product(z, u, Nx, Ny, DeltaT, DeltaX, DeltaY); // z = A * u
    vector_substract(r, b, z, n); // r = b - z
    memcpy(p, r, Nx*Ny*sizeof(double)); // p = r
    beta = sqrt(vector_scalar(r, r, n)); // beta = ||r||

    while (beta > eps && k < kmax) {

        matvect_product(z, p, Nx, Ny, DeltaT, DeltaX, DeltaY); // z = Ap

        alpha = vector_scalar(r, r, n) / vector_scalar(z, p, n); // alpha = <r,r> / <z,p>

        vector_coef(alphap, p, alpha, n); // alphax = alpha * p
        vector_sum(uplus, u, alphap, n); // x = x + alphap

        vector_coef(alphaz, z, alpha, n); // alphar = alpha * z
        vector_substract(rplus, r, alphaz, n); // r = r - alphaz

        gamma = vector_scalar(rplus, rplus, n) / vector_scalar(r, r, n); // gamma = <rplus,rplus> / <r,r>

        vector_coef(gammap, p, gamma, n); // gammap = gamma*p
        vector_sum(p, rplus, gammap, n); // p = r + gammap

        beta = sqrt(vector_scalar(rplus, rplus, n)); // beta = ||r||

        memcpy(r, rplus, Nx*Ny*sizeof(double));
        memcpy(u, uplus, Nx*Ny*sizeof(double));
        k++;
    }
    if (k > kmax) {
        printf("Tolérence non atteinte norme=%f", beta);
    }
    free(r), free(p), free(z), free(uplus), free(rplus), free(alphap), free(alphaz), free(gammap);
}
