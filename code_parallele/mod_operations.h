#include "/home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/parameter_and_display/mod_parameter.h"

void vector_sum(double* x, double* u, double* v, int n);

void vector_substract(double* x, double* u, double* v, int n);

void vector_coef(double* x, double* u, double c, int n);

double vector_scalar(double* u, double* v, struct data dt);

void matvect_product(double* b, double* x, struct data dt);

double Calculate_error(double* u_exact, double* u, struct data dt);