#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mod_operations.h"

double* gradient_conjugate(double* A, double* x, double* b, double eps, double kmax, int n) 
{
    int k = 0;

    double alpha, gamma, beta;
    double* r = (double*)malloc(n*sizeof(double));
    double* p = (double*)malloc(n*sizeof(double));
    double* z = (double*)malloc(n*sizeof(double));
    double* xplus = (double*)malloc(n*sizeof(double));
    double* rplus = (double*)malloc(n*sizeof(double));

    r = vector_substract(b, matvect_product(A, x, n), n);
    p = r;
    beta = sqrt(vector_scalar(r, r, n));
    
    printf("Before algo x = (");
    for (int i=0; i<n-1; i++) {
        printf("%f, ", x[i]);
    }
    printf("%f)\n", x[n-1]);
    printf("Norme = %f\n\n", beta);

    while (beta > eps && k < kmax) {
        z = matvect_product(A, p, n);
        alpha = vector_scalar(r, r, n) / vector_scalar(z, p, n);
        xplus = vector_sum(x, vector_coef(p, alpha, n), n);
        rplus = vector_substract(r, vector_coef(z, alpha, n), n);
        gamma = vector_scalar(rplus, rplus, n) / vector_scalar(r, r, n);
        p = vector_sum(rplus, vector_coef(p, gamma, n), n);
        beta = sqrt(vector_scalar(rplus, rplus, n));
        r = rplus;
        x = xplus;
        k += 1;

        printf("While algo x = (");
        for (int i=0; i<n-1; i++) {
            printf("%f, ", x[i]);
        }
        printf("%f)\n", x[n-1]);

        printf("Norme = %f , Itérations = %d \n\n", beta, k);
    }
    if (k > kmax) {
        printf("Tolérence non atteinte norme=%f", beta);
    }

    return x;
}