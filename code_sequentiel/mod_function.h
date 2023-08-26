#include "/home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/parameter_and_display/mod_parameter.h"

double func(double x, double y, double t, struct data dt);

double g(double x, double y, struct data dt);

double h(double x, double y, struct data dt);

double h_left(double y, double t, struct data dt);

double g_right(double y, double t, struct data dt);

double h_top(double x, double t, struct data dt);

double g_bottom(double x, double t, struct data dt);

double func_exact(double x, double y, double t, struct data dt);

double sol_init(double x, double y, struct data dt);