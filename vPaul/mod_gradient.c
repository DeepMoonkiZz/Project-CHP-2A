#include "mod_gradient.h"
#include "mod_matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

double norme(double *vect, int Imax){

  double val_norme = 0.0;

  for (int i = 0; i < Imax; i++){
    val_norme += vect[i]*vect[i];
  }

  return sqrt(val_norme);
}

double scalar(double *u, double *v, int Imax){

  double val_scalar = 0.0;

  for (int i = 0; i < Imax; i++){
    val_scalar += u[i+4]*v[i];
  }

  return val_scalar;
}

void conjugate_gradiant(double* u_1, double* u_0, int kmax, double D, double t, double dt, double dx, double dy, double eps, int Nx, int Ny){

    // NOTES :
    // u_1 désigne le vecteur u(n+1) donc celui à trouver (x dans Ax = b)
    // u_0 désigne le vecteur u(n)
    // L'équation est donc : A*u_(n+1) = u(n)
    //
    // Pour nos notations :
    // u(n) permet de remplir le vecteur b
    //
    // Les vecteurs concernés par la "marge de zéro" dans cette fonction void
    // sont : p, u_1 et u_0. On va donc les initialiser en amont pour ne pas le
    // faire à chaque itération.

    int Imax = Nx*Ny;
    double alpha, gamma, scalar_r;
    double * z = (double*)malloc(Imax*sizeof(double));
    double * r = (double*)malloc(Imax*sizeof(double));
    double * p = (double*)malloc((Imax+8)*sizeof(double));
    double * vect_b =(double*)malloc(Imax*sizeof(double));
    double * Au = (double*)malloc(Imax*sizeof(double));
    double new_p;

    build_vector_Au(dx,dy,u_0,Au,Nx,Ny);
    build_vector_b(dt,D,dx,dy,Nx,Ny,u_0,vect_b); //    r0 <- b - Ax0
                                                  //    p <- r0
    for (int i = 0; i < Imax+8; i++){
      p[i] = 0.0;
    }

    for (int i = 0; i < Imax; i++){
      r[i] = vect_b[i] - Au[i];
      p[i+4] = r[i];
    }

    free(vect_b);
    free(Au);

    double beta = norme(r,Imax);                //    Beta <- norme(r0)
    int k = 0;

    while((beta > eps) && (k <= kmax)){           //    Condition d'arrêt

      beta = norme(r,Imax);                     //    beta = norme(r)
      build_vector_Au(dx, dy, p, z, Nx, Ny);      //    z = Ap
      scalar_r = scalar(r,r,Imax);              //    <r,r>
      alpha = scalar_r/scalar(p,z,Imax);          //    alpha = <r,r>/<p,z>>


      for (int i = 0; i < Imax; i++){
        u_1[i+4] += alpha*p[i+4];                 //    u = u + alpha*p
      }

      for (int i = 0; i < Imax; i++){
        r[i] -= alpha*z[i];                     //    r = r - alpha*z
      }

      gamma = scalar(r,r,Imax)/scalar_r;

      for (int i = 0; i < Imax; i++){
        new_p = p[i+4];
        p[i+4] = r[i] + gamma*new_p;           // p = r + gamma*p
      }

      k++;
    }

    if (k > kmax){
      printf("Tolérance non atteinte: %lf\n", beta);
    }

    free(z);
    free(r);
    free(p);

}
