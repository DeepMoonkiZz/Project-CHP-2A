#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* vector_sum(double* u, double* v, int n)
// Sum of u+v = x
{
    double* x = (double*)malloc(n*sizeof(double));

    for (int i=0; i<n; i++) {
        x[i] = u[i] + v[i];
    }

    return x;
}


double* vector_substract(double* u, double* v, int n)
// Substraction of u-v = x
{
    double* x = (double*)malloc(n*sizeof(double));

    for (int i=0; i<n; i++) {
        x[i] = u[i] - v[i];
    }

    return x;
}


double* vector_coef(double* u, double c, int n)
// Apply coefficient to c*u = x
{
    double* x = (double*)malloc(n*sizeof(double));

    for (int i=0; i<n; i++) {
        x[i] = c * u[i];
    }

    return x;
}


double vector_scalar(double* u, double* v, int n)
// Dot product of u and v
{
    double x = 0;

    for (int i=0; i<n; i++) {
        x += u[i] * v[i];
    }

    return x;
}


double* matvect_product(double* A, double* x, int n) 
// Matrix vector product of A (matrix) and b (vector)
{
    double* b = (double*)malloc(n*sizeof(double));

    for (int i=0; i<n; i++) {
        b[i] = 0;
        for (int j=0; j<n; j++) {
            b[i] += A[i*n+j] * x[j];
        }
    }

    return b;
}