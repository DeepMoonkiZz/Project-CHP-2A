#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mod_gradient.h"
#include "mod_operations.h"


int main(int argc, char ** argv) 
{
    int n;

    double *A, *x, *b; 

    n = 10;

    A = (double*)malloc(n*n*sizeof(double));
    x = (double*)malloc(n*sizeof(double));
    b = (double*)malloc(n*sizeof(double));

    // Decla du vecteur x
    for (int i=0; i<n; i++) {
        x[i] = 0.;
        b[i] = i;
    }

    // Decla de la matrice M
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            A[i*n+j] = 0;
            if (i==j) {
                A[i*n+j] = 2;
            }
        }
    } 

    printf("Before algo x = (");
    for (int i=0; i<n-1; i++) {
        printf("%f, ", x[i]);
    }
    printf("%f)\n\n", x[n-1]);

    gradient_conjugate(A, x, b, 0.00001, 10000, n);

    printf("After algo x = (");
    for (int i=0; i<n-1; i++) {
        printf("%f, ", x[i]);
    }
    printf("%f)\n\n", x[n-1]);


    return 0;
}