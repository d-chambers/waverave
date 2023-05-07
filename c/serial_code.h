// serial_code.h
#ifndef SERIAL_CODE_H
#define SERIAL_CODE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void awe_2d_explicit_solver_homogeneous_8th_order(double **UUo, double **UUm, double dx, double dy, double dt, double v, double *F, int it, double sx, double sy, int nx, int ny);
int main();

#endif // SERIAL_CODE_H

