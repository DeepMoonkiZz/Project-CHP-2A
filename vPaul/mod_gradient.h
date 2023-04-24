#ifndef MOD_GRADIENT_H
#define MOD_GRADIENT_H

void conjugate_gradiant(double* u_1, double* u_0, int kmax, double D, double t, double dt, double dx, double dy, double eps, int Nx, int Ny);

double scalar(double *u, double *v, int Imax);

double norme(double *vect, int Imax);

#endif /* MOD_GRADIENT_H */
