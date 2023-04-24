#ifndef MOD_MATRIX_H
#define MOD_MATRIX_H

void build_vector_Au(double dx, double dy, double* u, double*v, int Nx, int Ny);

void build_vector_b(double dt, double D, double dx, double dy, int Nx, int Ny, double* u, double* vect_b);

void build_vector_Aub(double *u_n, double *u, double* Aub, double D, double dt, double dx, double dy, int Nx, int Ny);

#endif /* MOD_MATRIX_H */
