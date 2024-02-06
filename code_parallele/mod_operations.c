#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <openmpi/mpi.h>

#include "mod_operations.h"


void vector_sum(double* x, double* u, double* v, int n)
// Sum of u+v = x
{
    for (int i = 0; i < n; i++) {
        x[i] = u[i] + v[i];
    }
}


void vector_substract(double* x, double* u, double* v, int n)
// Substraction of u-v = x
{
    for (int i = 0; i < n; i++) {
        x[i] = u[i] - v[i];
    }
}


void vector_coef(double* x, double* u, double c, int n)
// Apply coefficient to c*u = x
{
    for (int i = 0; i < n; i++) {
        x[i] = c * u[i];
    }
}


double vector_scalar(double* u, double* v, struct data dt)
// Dot product of u and v
{
    double x = 0;
    for (int i = 0; i < dt.npart; i++) {
        x += u[i] * v[i];
    }
    MPI_Send(&x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    if (dt.rank == 0) {
        double val = 0, sum = 0;
        for (int i=0; i<dt.nproc; i++) {
            MPI_Recv(&val, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += val;
        }
        for (int i=0; i<dt.nproc; i++) {
            MPI_Send(&sum, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
    }
    MPI_Recv(&x, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return x;
}


void matvect_product(double* b, double* x, struct data dt)
// Matrix vector product of A (laplacian matrix) and x (vector) giving the result b
{
    double C, Cx, Cy;

    Cx = -dt.DeltaT*dt.D/(dt.DeltaX*dt.DeltaX);
    Cy = -dt.DeltaT*dt.D/(dt.DeltaY*dt.DeltaY);
    C = 1 - 2*Cx - 2*Cy;

    double ipx, imx, ipy, imy;

    if (dt.rank != 0) {
        // Send right flux
        if (dt.iBeg%dt.Nx != 0) {
            MPI_Send(&x[0], 1, MPI_DOUBLE, dt.rank-1, dt.Nx*dt.Ny + dt.iBeg, MPI_COMM_WORLD);
        }
        // Send up flux
        for (int i = 0; i<dt.Nx; i++) {
            MPI_Send(&x[i], 1, MPI_DOUBLE, dt.rank-1, dt.iBeg + i, MPI_COMM_WORLD);
        }
    }
    if (dt.rank != dt.nproc - 1) {
        // Send left flux
        if ((dt.npart-1)%dt.Nx != dt.Nx-1) {
            MPI_Send(&x[dt.npart-1], 1, MPI_DOUBLE, dt.rank+1, dt.Nx*dt.Ny + dt.iEnd, MPI_COMM_WORLD);
        }
        // Send down flux
        for (int i = dt.npart-dt.Nx; i<dt.npart; i++) {
            MPI_Send(&x[i], 1, MPI_DOUBLE, dt.rank+1, dt.iBeg + i, MPI_COMM_WORLD);
        }
    }


    for (int i = dt.iBeg; i <= dt.iEnd; i++) {
        b[i-dt.iBeg] = C * x[i-dt.iBeg];

        // Right flux
        if (i%dt.Nx != dt.Nx-1) {
            if (i+1 > dt.iEnd && dt.rank != dt.nproc-1) {
                MPI_Recv(&ipx, 1, MPI_DOUBLE, dt.rank+1, dt.Nx*dt.Ny + i+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
                ipx = x[i-dt.iBeg+1];
            }
            b[i-dt.iBeg] += Cx * ipx;
        }
        // Left flux
        if (i%dt.Nx != 0) {        
            if (i-1 < dt.iBeg && dt.rank != 0) {
                MPI_Recv(&imx, 1, MPI_DOUBLE, dt.rank-1, dt.Nx*dt.Ny + i-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
                imx = x[i-dt.iBeg-1];
            }
            b[i-dt.iBeg] += Cx * imx;
        }
        // Up flux
        if (i/dt.Nx != dt.Ny-1) {
            if (i+dt.Nx > dt.iEnd && dt.rank != dt.nproc-1) {
                MPI_Recv(&ipy, 1, MPI_DOUBLE, dt.rank+1, i+dt.Nx, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
                ipy = x[i-dt.iBeg+dt.Nx];
            }
            b[i-dt.iBeg] += Cy * ipy;
        }
        // Down flux
        if (i/dt.Nx != 0) {
            if (i-dt.Nx < dt.iBeg && dt.rank != 0) {
                MPI_Recv(&imy, 1, MPI_DOUBLE, dt.rank-1, i-dt.Nx, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
                imy = x[i-dt.iBeg-dt.Nx];
            }
            b[i-dt.iBeg] += Cy * imy;
        }
    }

    if (dt.function==4 || dt.function==5) {
        for (int i = dt.iBeg; i <= dt.iEnd; i++) {
            if (i/dt.Nx == dt.Ny-1) {
                b[i-dt.iBeg] += Cy * x[i-dt.iBeg];
            }
            if (i%dt.Nx == 0) {
                b[i-dt.iBeg] += Cx * x[i-dt.iBeg]; 
            }
        }
    }
}


double Calculate_error(double* u_exact, double* u, struct data dt)
{
    double sum=0;
    for (int i=0; i<dt.npart; i++) {
        sum += (u_exact[i]-u[i])*(u_exact[i]-u[i]);
    }
    MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);    
    if (dt.rank == 0) {
        double val = 0;
        sum = 0;
        for (int i=0; i<dt.nproc; i++) {
            MPI_Recv(&val, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += val;
        }
    }
    return sum/(dt.Nx*dt.Ny);
}


/*
void matvect_product(double* b, double* x, struct data dt)
// Matrix vector product of A (laplacian matrix) and x (vector) giving the result b
{
    double C, Cx, Cy;

    Cx = -dt.DeltaT*dt.D/(dt.DeltaX*dt.DeltaX);
    Cy = -dt.DeltaT*dt.D/(dt.DeltaY*dt.DeltaY);
    C = 1 - 2*Cx - 2*Cy;

    double ipx, imx, ipy, imy;

    MPI_Request send_request;
    MPI_Request recv_request;

    if (dt.rank != 0) {
        // Send right flux
        if (dt.iBeg%dt.Nx != 0) {
            MPI_Isend(&x[0], 1, MPI_DOUBLE, dt.rank-1, dt.Nx*dt.Ny + dt.iBeg, MPI_COMM_WORLD, &send_request);
        }
        // Send up flux
        for (int i = 0; i<dt.Nx; i++) {
            MPI_Isend(&x[i], 1, MPI_DOUBLE, dt.rank-1, dt.iBeg + i, MPI_COMM_WORLD, &send_request);
        }
    }
    if (dt.rank != dt.nproc - 1) {
        // Send left flux
        if ((dt.npart-1)%dt.Nx != dt.Nx-1) {
            MPI_Isend(&x[dt.npart-1], 1, MPI_DOUBLE, dt.rank+1, dt.Nx*dt.Ny + dt.iEnd, MPI_COMM_WORLD, &send_request);
        }
        // Send down flux
        for (int i = dt.npart-dt.Nx; i<dt.npart; i++) {
            MPI_Isend(&x[i], 1, MPI_DOUBLE, dt.rank+1, dt.iBeg + i, MPI_COMM_WORLD, &send_request);
        }
    }
        MPI_Barrier(MPI_COMM_WORLD);


    for (int i = dt.iBeg; i <= dt.iEnd; i++) {
        b[i-dt.iBeg] = C * x[i-dt.iBeg];

        // Right flux
        if (i%dt.Nx != dt.Nx-1) {
            if (i+1 > dt.iEnd && dt.rank != dt.nproc-1) {
                MPI_Irecv(&ipx, 1, MPI_DOUBLE, dt.rank+1, dt.Nx*dt.Ny + i+1, MPI_COMM_WORLD, &recv_request);
                MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
            }
            else {
                ipx = x[i-dt.iBeg+1];
            }
            b[i-dt.iBeg] += Cx * ipx;
        }
        // Left flux
        if (i%dt.Nx != 0) {        
            if (i-1 < dt.iBeg && dt.rank != 0) {
                MPI_Irecv(&imx, 1, MPI_DOUBLE, dt.rank-1, dt.Nx*dt.Ny + i-1, MPI_COMM_WORLD, &recv_request);
                MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
            }
            else {
                imx = x[i-dt.iBeg-1];
            }
            b[i-dt.iBeg] += Cx * imx;
        }
        // Up flux
        if (i/dt.Nx != dt.Ny-1) {
            if (i+dt.Nx > dt.iEnd && dt.rank != dt.nproc-1) {
                MPI_Irecv(&ipy, 1, MPI_DOUBLE, dt.rank+1, i+dt.Nx, MPI_COMM_WORLD, &recv_request);
                MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
            }
            else {
                ipy = x[i-dt.iBeg+dt.Nx];
            }
            b[i-dt.iBeg] += Cy * ipy;
        }
        // Down flux
        if (i/dt.Nx != 0) {
            if (i-dt.Nx < dt.iBeg && dt.rank != 0) {
                MPI_Irecv(&imy, 1, MPI_DOUBLE, dt.rank-1, i-dt.Nx, MPI_COMM_WORLD, &recv_request);
                MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
            }
            else {
                imy = x[i-dt.iBeg-dt.Nx];
            }
            b[i-dt.iBeg] += Cy * imy;
        }
    }

    if (dt.function==4 || dt.function==5) {
        for (int i = dt.iBeg; i <= dt.iEnd; i++) {
            if (i/dt.Nx == dt.Ny-1) {
                b[i-dt.iBeg] += Cy * x[i-dt.iBeg];
            }
            if (i%dt.Nx == 0) {
                b[i-dt.iBeg] += Cx * x[i-dt.iBeg]; 
            }
        }
    }
}
*/




/*
   MPI_Request send_request[2];
    MPI_Request recv_requests[2];

    double ipx, imx, ipy, imy;
    double *proc_plus, *proc_moins;
    proc_plus = (double*)malloc(dt.Nx*sizeof(double));
    proc_moins = (double*)malloc(dt.Nx*sizeof(double));
    printf("\n\n");
                for (int i=0; i < dt.npart; i++) {
                printf("x i: %i procm: %f\n", i, x[i]);
            }

    printf("\n\n");
    if (dt.nproc != 1) {
        if (dt.rank == 0) {
            for (int i=0; i < dt.Nx; i++) {
                proc_moins[i] = x[dt.npart - dt.Nx + i];
                //printf("i: %i procm: %f\n", i, x[dt.npart - dt.Nx + i]);
            }
            MPI_Isend(proc_moins, dt.Nx, MPI_DOUBLE, dt.rank + 1, dt.rank, MPI_COMM_WORLD, &send_request[0]);
            MPI_Irecv(proc_plus, dt.Nx, MPI_DOUBLE, dt.rank + 1, dt.nproc*(dt.rank + 1), MPI_COMM_WORLD, &recv_requests[1]);
            MPI_Wait(&recv_requests[1], MPI_STATUS_IGNORE);
        }
        else if (dt.rank == dt.nproc - 1) {
            for (int i=0; i < dt.Nx; i++) {
                proc_plus[i] = x[i];
            }
            MPI_Isend(proc_plus, dt.Nx, MPI_DOUBLE, dt.rank - 1, dt.rank*dt.nproc, MPI_COMM_WORLD, &send_request[1]);
            MPI_Irecv(proc_moins, dt.Nx, MPI_DOUBLE, dt.rank - 1, dt.rank - 1, MPI_COMM_WORLD, &recv_requests[0]);
            MPI_Wait(&recv_requests[0], MPI_STATUS_IGNORE);
        }
        else {
            for (int i=0; i < dt.Nx; i++) {
                proc_plus[i] = x[i];
                proc_moins[i] = x[dt.npart - dt.Nx + i];
            }
            MPI_Isend(proc_moins, dt.Nx, MPI_DOUBLE, dt.rank + 1, dt.rank, MPI_COMM_WORLD, &send_request[0]);
            MPI_Isend(proc_plus, dt.Nx, MPI_DOUBLE, dt.rank - 1, dt.rank*dt.nproc, MPI_COMM_WORLD, &send_request[1]);
            free(proc_plus), free(proc_moins);

            MPI_Irecv(proc_moins, dt.Nx, MPI_DOUBLE, dt.rank - 1, dt.rank - 1, MPI_COMM_WORLD, &recv_requests[1]);
            MPI_Irecv(proc_plus, dt.Nx, MPI_DOUBLE, dt.rank + 1, dt.nproc*(dt.rank + 1), MPI_COMM_WORLD, &recv_requests[0]);

            MPI_Wait(&recv_requests[0], MPI_STATUS_IGNORE);
            MPI_Wait(&recv_requests[1], MPI_STATUS_IGNORE);
        }
    }

    for (int i = dt.iBeg; i <= dt.iEnd; i++) {
        b[i-dt.iBeg] = C * x[i-dt.iBeg];

        // Right flux
        if (i%dt.Nx != dt.Nx-1) {
            if (i+1 > dt.iEnd && i%dt.Nx+1 < dt.Nx) {
                ipx = proc_plus[0];
            }
            else {
                ipx = x[i-dt.iBeg+1];
            }
            b[i-dt.iBeg] += Cx * ipx;
        }
        // Left flux
        if (i%dt.Nx != 0) {        
            if (i-1 < dt.iBeg && i%dt.Nx-1 >= 0) {
                imx = proc_moins[dt.Nx - 1];
            }
            else {
                imx = x[i-dt.iBeg-1];
            }
            b[i-dt.iBeg] += Cx * imx;
        }
        // Up flux
        if (i/dt.Nx != dt.Ny-1) {
            if (i+dt.Nx > dt.iEnd) {
                ipy = proc_plus[(dt.iEnd+i)%dt.Nx];
            }
            else {
                ipy = x[i-dt.iBeg+dt.Nx];
            }
            b[i-dt.iBeg] += Cy * ipy;
        }
        // Down flux
        if (i/dt.Nx != 0) {
            if (i-dt.Nx < dt.iBeg && (i-dt.Nx)/dt.Nx >= 0) {
                imy = proc_moins[i-dt.iBeg];
            }
            else {
                imy = x[i-dt.iBeg-dt.Nx];
            }
            b[i-dt.iBeg] += Cy * imy;
        }
    }

    if (dt.function==4 || dt.function==5) {
        for (int i = dt.iBeg; i <= dt.iEnd; i++) {
            if (i/dt.Nx == dt.Ny-1) {
                b[i-dt.iBeg] += Cy * x[i-dt.iBeg];
            }
            if (i%dt.Nx == 0) {
                b[i-dt.iBeg] += Cx * x[i-dt.iBeg]; 
            }
        }
    }
    free(proc_plus), free(proc_moins);
*/