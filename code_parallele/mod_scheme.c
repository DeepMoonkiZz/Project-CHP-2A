#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <openmpi/mpi.h>
#include "/home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/fonction/charge.h"

#include "mod_function.h"
#include "mod_scheme.h"


void Build_vect_b(double* b, double* u, double t, struct data dt)
{
    double Cx = dt.DeltaT*dt.D / (dt.DeltaX*dt.DeltaX);
    double Cy = dt.DeltaT*dt.D / (dt.DeltaY*dt.DeltaY);


    // Remplissage du vecteur b avec les composantes du vecteur U

    memcpy(b, u, dt.npart*sizeof(double));

    // Remplissage du vecteur b avec les composantes du vecteur F

    for (int i = dt.iBeg; i <= dt.iEnd; i++) {
        b[i-dt.iBeg] += dt.DeltaT * func(dt.xmin + (i%dt.Nx)*dt.DeltaX, dt.ymin + (i/dt.Nx)*dt.DeltaY, t, dt);
    }

    // Remplissage du vecteur b avec les composantes des conditions limites
    if (dt.function==4 || dt.function==5) {
        // Up and down side
        for (int i = 0; i < dt.Nx; i++) {
            if (i >= dt.iBeg && i <= dt.iEnd) {
                b[i-dt.iBeg] += Cy * g_bottom(dt.xmin + (i%dt.Nx)*dt.DeltaX, t-dt.DeltaT, dt);
            }
            if (dt.Nx*(dt.Ny-1) + i >= dt.iBeg && dt.Nx*(dt.Ny-1) + i <= dt.iEnd) {
                b[dt.Nx*(dt.Ny-1) + i - dt.iBeg] -= dt.DeltaY * Cy * h_top(dt.xmin + (i%dt.Nx)*dt.DeltaX, t-dt.DeltaT, dt);
            }
        }
        // Left and right side
        for (int j = 0; j < dt.Ny; j++) {
            if (dt.Nx*j >= dt.iBeg && dt.Nx*j <= dt.iEnd) {
                b[dt.Nx*j-dt.iBeg] -= dt.DeltaX * Cx * h_left(dt.ymin + j*dt.DeltaY, t-dt.DeltaT, dt); 
            }
            if (dt.Nx*j + dt.Nx-1 >= dt.iBeg && dt.Nx*j + dt.Nx-1 <= dt.iEnd) {
                b[dt.Nx*j + dt.Nx-1 - dt.iBeg] += Cx * g_right(dt.ymin + j*dt.DeltaY, t-dt.DeltaT, dt);
            }
        }
    }
    else {    
        // Up and down side    
        for (int i = 0; i < dt.Nx; i++) {
            if (i >= dt.iBeg && i <= dt.iEnd) {
                b[i-dt.iBeg] += Cy * g(dt.xmin + i*dt.DeltaX, dt.ymin - dt.DeltaY, dt);
            }
            if (dt.Nx*(dt.Ny-1) + i >= dt.iBeg && dt.Nx*(dt.Ny-1) + i <= dt.iEnd) {
                b[dt.Nx*(dt.Ny-1) + i - dt.iBeg] += Cy * g(dt.xmin + i*dt.DeltaX, dt.ymax + dt.DeltaY, dt);
            }
        }
        // Left and right side
        for (int j = 0; j < dt.Ny; j++) {
            if (dt.Nx*j >= dt.iBeg && dt.Nx*j <= dt.iEnd) {
                b[dt.Nx*j-dt.iBeg] += Cx * h(dt.xmin - dt.DeltaX, dt.ymin + j*dt.DeltaY, dt);
            }
            if (dt.Nx*j + dt.Nx-1 >= dt.iBeg && dt.Nx*j + dt.Nx-1 <= dt.iEnd) {
                b[dt.Nx*j + dt.Nx-1 - dt.iBeg] += Cx * h(dt.xmax + dt.DeltaX, dt.ymin + j*dt.DeltaY, dt);
            }
        }
    }
}


void Build_u_exact(double* u_exact, double t, struct data dt)
{
    for (int i = dt.iBeg; i <= dt.iEnd; i++) {
        u_exact[i-dt.iBeg] = func_exact(dt.xmin + (i%dt.Nx)*dt.DeltaX, dt.ymin + (i/dt.Nx)*dt.DeltaY, t, dt);
    }
}


void Build_comp_vect(double *u, double *u_part, struct data dt) 
{
    MPI_Request send_request;
    MPI_Request *recv_request;

    int iBeg, iEnd, n_part;
    double *y_part;

    // Send u_part and u_exact_part to proc 0
    MPI_Isend(u_part, dt.npart, MPI_DOUBLE, 0, dt.rank, MPI_COMM_WORLD, &send_request);

    if (dt.rank == 0) {
        // Receive u_part and build u
        recv_request = (MPI_Request*)malloc(sizeof(MPI_Request) * (dt.nproc));
        for (int i=0; i<dt.nproc; i++) {
            charge(i, dt.Nx*dt.Ny, dt.nproc, &iBeg, &iEnd);
            n_part = iEnd-iBeg+1;
            y_part = (double*)malloc(n_part*sizeof(double));
            
            MPI_Irecv(y_part, n_part, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &recv_request[i]);

            MPI_Wait(&recv_request[i], MPI_STATUS_IGNORE);
            
            for (int j=iBeg; j<=iEnd; j++) {
                u[j] = y_part[j-iBeg];
            }
            free(y_part);
        }
        free(recv_request);
    }
}