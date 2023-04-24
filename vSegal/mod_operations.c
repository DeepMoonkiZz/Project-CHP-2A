#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void vector_sum(double* u, double* v, double* x, int n)
// Sum of u+v = x
{
    for (int i=0; i<n; i++) {
        x[i] = u[i] + v[i];
    }
}


void vector_substract(double* u, double* v, double* x, int n)
// Substraction of u-v = x
{
    for (int i=0; i<n; i++) {
        x[i] = u[i] - v[i];
    }
}


void vector_coef(double* u, double c, double* x, int n)
// Apply coefficient to c*u = x
{
    for (int i=0; i<n; i++) {
        x[i] = c * u[i];
    }
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


void matvect_product(double* A, double* x, double* b, int n) 
// Matrix vector product of A (matrix) and b (vector)
{
    for (int i=0; i<n; i++) {
        b[i] = 0;
        for (int j=0; j<n; j++) {
            b[i] += A[i*n+j] * x[j];
        }
    }
}