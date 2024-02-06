#ifndef MOD_SCHEME_H
#define MOD_SCHEME_H

#include "/home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/parameter_and_display/mod_parameter.h"

void Build_vect_b(double* b, double* u, double t, struct data dt);

void Build_u_exact(double* u_exact, double t, struct data dt);

#endif // MOD_SCHEME_H