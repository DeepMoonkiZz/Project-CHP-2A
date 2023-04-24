#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mod_operations.h"

void gradient_conjugate(double* A, double* x, double* b, double eps, double kmax, int n) 
{
    int k = 0;

    double alpha, gamma, beta;
    double* r = (double*)malloc(n*sizeof(double));
    double* p = (double*)malloc(n*sizeof(double));
    double* z = (double*)malloc(n*sizeof(double));
    double* Ap = (double*)malloc(n*sizeof(double));
    double* alphap = (double*)malloc(n*sizeof(double));
    double* alphaz = (double*)malloc(n*sizeof(double));
    double* gammap = (double*)malloc(n*sizeof(double));

    matvect_product(A, x, Ap, n);
    vector_substract(b, Ap, r, n);
    p = r;
    beta = sqrt(vector_scalar(r, r, n));
    printf("Norme: %f\n", beta);

    while (beta > eps && k < kmax) {
        z = matvect_product(A, p, n);
        alpha = vector_scalar(r, r, n) / vector_scalar(z, p, n);
        x = vector_sum(x, vector_coef(p, alpha, n), n);
        gamma = vector_scalar(vector_substract(r, vector_coef(z, alpha, n), n), vector_substract(r, vector_coef(z, alpha, n), n), n) / vector_scalar(r, r, n);
        p = vector_sum(vector_substract(r, vector_coef(z, alpha, n), n), vector_coef(p, gamma, n), n);
        beta = sqrt(vector_scalar(r, r, n));
        r = vector_substract(r, vector_coef(z, alpha, n), n);
        k += 1;

        printf("While algo x = (");
        for (int i=0; i<n-1; i++) {
            printf("%f, ", x[i]);
        }
        printf("%f)\n\n", x[n-1]);

        printf("Norme = %f , Itérations = %d \n", beta, k);
    }
    if (k > kmax) {
        printf("Tolérence non atteinte norme=%f", beta);
    }
}