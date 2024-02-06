#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <openmpi/mpi.h>

#include "mod_operations.h"
#include "mod_gradient.h"

void gradient_conjugate(double *u, double *b, struct data dt)
{
    int k = 0;

    double alpha, gamma, beta;
    double *r, *p, *z, *uplus, *rplus, *alphap, *alphaz, *gammap;
    r = (double*)malloc(dt.npart*sizeof(double));
    p = (double*)malloc(dt.npart*sizeof(double));
    z = (double*)malloc(dt.npart*sizeof(double));
    uplus = (double*)malloc(dt.npart*sizeof(double));
    rplus = (double*)malloc(dt.npart*sizeof(double));
    alphap = (double*)malloc(dt.npart*sizeof(double));
    alphaz = (double*)malloc(dt.npart*sizeof(double));
    gammap = (double*)malloc(dt.npart*sizeof(double));

    matvect_product(z, u, dt); // z = A * u
        
    vector_substract(r, b, z, dt.npart); // r = b - z
    memcpy(p, r, dt.npart*sizeof(double)); // p = r
    beta = sqrt(vector_scalar(r, r, dt)); // beta = ||r||

    while (beta > dt.eps && k < dt.Kmax) {

        matvect_product(z, p, dt); // z = Ap
        alpha = vector_scalar(r, r, dt) / vector_scalar(z, p, dt); // alpha = <r,r> / <z,p>

        vector_coef(alphap, p, alpha, dt.npart); // alphax = alpha * p
        vector_sum(uplus, u, alphap, dt.npart); // x = x + alphap

        vector_coef(alphaz, z, alpha, dt.npart); // alphar = alpha * z
        vector_substract(rplus, r, alphaz, dt.npart); // r = r - alphaz

        gamma = vector_scalar(rplus, rplus, dt) / vector_scalar(r, r, dt); // gamma = <rplus,rplus> / <r,r>

        vector_coef(gammap, p, gamma, dt.npart); // gammap = gamma*p
        vector_sum(p, rplus, gammap, dt.npart); // p = r + gammap

        beta = sqrt(vector_scalar(rplus, rplus, dt)); // beta = ||r||

        memcpy(r, rplus, dt.npart*sizeof(double));
        memcpy(u, uplus, dt.npart*sizeof(double));
        k++;
    }
    if (k > dt.Kmax) {
        printf("Tolérence non atteinte norme=%f", beta);
    }
    free(r);
    free(p);
    free(z);
    free(uplus);
    free(rplus);
    free(alphap);
    free(alphaz); 
    free(gammap);
}