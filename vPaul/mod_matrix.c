#include "mod_matrix.h"
#include "mod_function.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void build_vector_Au(double dx, double dy, double* u, double*v, int Nx, int Ny){

  // NOTES :
  // u désigne le vecteur dans Au
  // v désigne le résultat v = A*u

  // MÉTHODE DE CACLCUL : pour éviter les problèmes de bords dans la matrice,
  // on va travailler avec des vecteurs u à Imax+8 cases remplies de zéros.
  // Cela permettra de créer une "marge de zéros" pour la multiplication.
  // On va donc décaler les indices de +4 dans CHAQUE calcul comprenant un des
  // vecteurs u concernés. :)

  double b = -1./(dx*dx);
  double c = -1./(dy*dy);
  double a = 1. + 2*b + 2*c;
  int Imax = Nx*Ny;

  for (int k = 0; k < Imax; k++){

      if((k+1)%(Nx) != 0){

          v[k] = c*u[k] + b*u[k+3] + a*u[k+4] + b*u[k+5] + c*u[k+8];
      }
      else{

          v[k] = c*u[k] + b*u[k+3] + a*u[k+4] + c*u[k+8];
          k++;
          v[k] = c*u[k] + a*u[k+4] + b*u[k+5] + c*u[k+8];
      }
  }
}

void build_vector_b(double dt, double D, double dx, double dy, int Nx, int Ny, double* u, double* vect_b){

  // NOTES :
  // u désigne le vecteur u_n
  // vect_b désigne le second membre b = u_n + F_n

  int Imax = Nx*Ny;
  int k;
  double b = -1./(dx*dx);
  double c = -1./(dy*dy);
  double * F = malloc(Nx*Ny*sizeof(double));
  double * BC = malloc(Nx*Ny*sizeof(double));
  double * G = malloc(Nx*Ny*sizeof(double));

  // On remplit d'abord F_n
  // La fonction peut changer, attention aux arguments de build_vector_b

  for (int j = 1; j < Ny+1; j++){
    for (int i = 1; i < Nx+1; i++){
      k = Nx*(j-1) + (i-1);
      F[k] = f_1(i*dx,j*dy);
    }
  }

  // On va maintenant remplir BC_n
  // Les fonctions h et g peuvent aussi changer

  // On traite d'abord les 4 coins :
  BC[0] = 0.0;
  BC[Nx-1] = 0.0;
  BC[Nx*(Ny-1)] = 0.0;
  BC[Nx*(Ny-1)+(Nx-1)] = 0.0;

  // On traite maintenant les bords Sigma_0 :
  for (int i = 2; i < Nx; i++){
    BC[i-1] = 0.0; // j = 1
    BC[Nx*(Ny-1)+(i-1)] = 0.0; // j = Ny
  }

  // Enfin, on traite les bords Sigma_1
  for (int j = 2; j < Ny; j++){
    BC[Nx*(j-1)] = 0.0; // j = 1
    BC[Nx*(Ny-1)+(Nx-1)] = 0.0; // j = Ny
  }

  // On peut maintenant calculer G_n = F_n - BC_n
  for (int i = 0; i < Imax; i++){
    G[i] = F[i] - BC[i];
  }

  free(BC);
  free(F);

  // Enfin, on calcul b = 1/dt*u + 1/D*G
  for (int i = 0; i < Imax; i ++){
    vect_b[i] = (1./dt)*u[i] + (1./D)*G[i];
  }

  free(G);

}

void build_vector_Aub(double *u_n, double *u, double* Aub, double D, double dt, double dx, double dy, int Nx, int Ny){

  // NOTES :
  // u_n désigne le vecteur u_n contenu dans b
  // u désigne le vecteur de Au
  // Aub désigne le résultat de Aub = Au - b

  int Imax = Nx*Ny;

  // Construction du vecteur b
  double * vect_b = malloc(Imax*sizeof(double));
  build_vector_b(dt,D,dx,dy,Nx,Ny,u_n,vect_b);

  // construction du vecteur Au
  double * Au = malloc(Imax*sizeof(double));
  build_vector_Au(dx,dy,u,Au,Nx,Ny);

  // Simple différence pour calculer Au - b
  for(int i = 0; i < Imax; i++){

    Aub[i+4] = vect_b[i] - Au[i];
  }

  free(Au);
  free(vect_b);

}
