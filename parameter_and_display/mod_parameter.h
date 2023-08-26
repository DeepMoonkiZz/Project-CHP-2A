#ifndef MOD_PARAMETER_H
#define MOD_PARAMETER_H

#include <stdio.h>

struct data
{
    int function;
    int Kmax;
    double eps;
    double Tmax;

    double xmax;
    double xmin;
    double ymax; 
    double ymin;

    int Nx;
    int Ny;
    double Lx;
    double Ly;

    double D;

    double DeltaT;
    double DeltaX;
    double DeltaY;

    int rank;
    int nproc;    
    int iBeg;
    int iEnd;
};

void Read_parameter(struct data *dt);

#endif // MOD_PARAMETER_H