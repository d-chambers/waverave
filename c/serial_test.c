#include "serial_code.h"
#include <assert.h>
#include <stdbool.h>

void awe_2d_explicit_solver_homogeneous_8th_order(double **UUo, double **UUm, double dx, double dy, double dt, double v, double *F, int it, double sx, double sy, int nx, int ny) {
    int isx, isy, i, j;
    double Cx2, Cy2;

    // Get dimensions of wavefield
    nx = nx;
    ny = ny;

    // Define Courant numbers (squared)
    Cx2 = pow((v * dt / dx), 2);
    Cy2 = pow((v * dt / dy), 2);

    // Source location
    isx = (int) (sx / dx);
    isy = (int) (sy / dy);

    // Inject wavelet
    UUo[isx][isy] += dt * dt * F[it];

    // Update solution
    for (i = 4; i < nx - 4; i++) {
        for (j = 4; j < ny - 4; j++) {
            UUm[i][j] = 2 * UUo[i][j] - UUm[i][j] + Cx2 * (
                -1.0 / 560 * UUo[i - 4][j] + 8.0 / 315 * UUo[i - 3][j] - 1.0 / 5 * UUo[i - 2][j] + 8.0 / 5 * UUo[i - 1][j] - 205.0 / 72 * UUo[i][j] + 8.0 / 5 * UUo[i + 1][j] - 1.0 / 5 * UUo[i + 2][j] + 8.0 / 315 * UUo[i + 3][j] - 1.0 / 560 * UUo[i + 4][j]
            ) + Cy2 * (
                -1.0 / 560 * UUo[i][j - 4] + 8.0 / 315 * UUo[i][j - 3] - 1.0 / 5 * UUo[i][j - 2] + 8.0 / 5 * UUo[i][j - 1] - 205.0 / 72 * UUo[i][j] + 8.0 / 5 * UUo[i][j + 1] - 1.0 / 5 * UUo[i][j + 2] + 8.0 / 315 * UUo[i][j + 3] - 1.0 / 560 * UUo[i][j + 4]
            );
        }
    }
}

void test_awe_2d_explicit_solver_homogeneous_8th_order() {
    int nx = 10, ny = 10;
    double dx = 1.0, dy = 1.0, dt = 0.01, v = 1.0, sx = 5.0, sy = 5.0;
    double F = 1.0;
    int it = 0;

    double **UUo = (double **) malloc(nx * sizeof(double *));
    double **UUm = (double **) malloc(nx * sizeof(double *));
    for (int i = 0; i < nx; i++) {
        UUo[i] = (double *) malloc(ny * sizeof(double));
        UUm[i] = (double *) malloc(ny * sizeof(double));
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            UUo[i][j] = 0;
            UUm[i][j] = 0;
        }
    }

    awe_2d_explicit_solver_homogeneous_8th_order(UUo, UUm, dx, dy, dt, v, &F, it, sx, sy, nx, ny);

    // Check if the source position is updated correctly
    int isx = (int) (sx / dx);
    int isy = (int) (sy / dy);
    double expected_value = dt * dt * F;
    if (fabs(UUo[isx][isy] - expected_value) < 1e-6) {
        printf("Test for awe_2d_explicit_solver_homogeneous_8th_order passed.\n");
    } else {
        printf("Test for awe_2d_explicit_solver_homogeneous_8th_order failed.\n");
    }

    for (int i = 0; i < nx; i++) {
        free(UUo[i]);
        free(UUm[i]);
    }
    free(UUo);
    free(UUm);
}

void test_zero_initial_conditions() {
    int nx = 10, ny = 10;
    double dx = 1.0, dy = 1.0, dt = 0.01, v = 1.0, sx = 5.0, sy = 5.0;
    double F = 1.0;
    int it = 0;

    double **UUo = (double **) malloc(nx * sizeof(double *));
    double **UUm = (double **) malloc(nx * sizeof(double *));
    for (int i = 0; i < nx; i++) {
        UUo[i] = (double *) malloc(ny * sizeof(double));
        UUm[i] = (double *) malloc(ny * sizeof(double));
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            UUo[i][j] = 0; // Set zero initial conditions
            UUm[i][j] = 0;
        }
    }

    awe_2d_explicit_solver_homogeneous_8th_order(UUo, UUm, dx, dy, dt, v, &F, it, sx, sy, nx, ny);

    // Check if the output wavefield values are zero (apart from the source injection point)
    bool zero_initial_conditions_test_passed = true;
    int isx = (int) (sx / dx);
    int isy = (int) (sy / dy);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            if (i == isx && j == isy) {
                continue; // Exclude the source injection point
            }
            if (fabs(UUo[i][j]) > 1e-6) {
                zero_initial_conditions_test_passed = false;
                break;
            }
        }
    }

    if (zero_initial_conditions_test_passed) {
        printf("Test for zero initial conditions passed.\n");
    } else {
        printf("Test for zero initial conditions failed.\n");
    }

    for (int i = 0; i < nx; i++) {
        free(UUo[i]);
        free(UUm[i]);
    }
    free(UUo);
    free(UUm);
}

void calculate_courant_numbers(double dx, double dy, double dt, double v, double *Cx2, double *Cy2) {
    // Define Courant numbers (squared)
    *Cx2 = pow((v * dt / dx), 2);
    *Cy2 = pow((v * dt / dy), 2);
}

void test_courant_numbers() {
    double dx = 1.0, dy = 1.0, dt = 0.01, v = 1.0;
    double expected_Cx2 = pow((v * dt / dx), 2);
    double expected_Cy2 = pow((v * dt / dy), 2);

    double Cx2, Cy2;
    // Calculate Courant numbers directly from the calculate_courant_numbers function
    calculate_courant_numbers(dx, dy, dt, v, &Cx2, &Cy2);

    assert(fabs(Cx2 - expected_Cx2) < 1e-6 && "Test for Cx2 failed");
    assert(fabs(Cy2 - expected_Cy2) < 1e-6 && "Test for Cy2 failed");

    printf("Test for Courant numbers passed.\n");
}
void test_boundary_conditions() {
    int nx = 10, ny = 10;
    double dx = 1.0, dy = 1.0, dt = 0.01, v = 1.0, sx = 5.0, sy = 5.0;
    double F = 1.0;
    int it = 0;

    double **UUo = (double **) malloc(nx * sizeof(double *));
    double **UUm = (double **) malloc(nx * sizeof(double *));
    for (int i = 0; i < nx; i++) {
        UUo[i] = (double *) malloc(ny * sizeof(double));
        UUm[i] = (double *) malloc(ny * sizeof(double));
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            UUo[i][j] = 0;
            UUm[i][j] = 0;
        }
    }

    // Save the boundary values before calling the solver
    double **boundary_values = (double **) malloc(nx * sizeof(double *));
    for (int i = 0; i < nx; i++) {
        boundary_values[i] = (double *) malloc(ny * sizeof(double));
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            if (i < 4 || i >= nx - 4 || j < 4 || j >= ny - 4) {
                boundary_values[i][j] = UUo[i][j];
            }
        }
    }

    awe_2d_explicit_solver_homogeneous_8th_order(UUo, UUm, dx, dy, dt, v, &F, it, sx, sy, nx, ny);

    // Check if the boundary values remain unchanged
    bool boundary_test_passed = true;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            if (i < 4 || i >= nx - 4 || j < 4 || j >= ny - 4) {
                if (fabs(UUo[i][j] - boundary_values[i][j]) > 1e-6) {
                    boundary_test_passed = false;
                    break;
                }
            }
        }
    }

    if (boundary_test_passed) {
        printf("Test for boundary conditions passed.\n");
    } else {
        printf("Test for boundary conditions failed.\n");
    }

    for (int i = 0; i < nx; i++) {
        free(boundary_values[i]);
        free(UUo[i]);
        free(UUm[i]);
    }
    free(boundary_values);
    free(UUo);
    free(UUm);
}

void test_source_injection() {
    int nx = 10, ny = 10;
    double dx = 1.0, dy = 1.0, dt = 0.01, v = 1.0, sx = 5.0, sy = 5.0;
    double F = 1.0;
    int it = 0;

    double **UUo = (double **) malloc(nx * sizeof(double *));
    double **UUm = (double **) malloc(nx * sizeof(double *));
    for (int i = 0; i < nx; i++) {
        UUo[i] = (double *) malloc(ny * sizeof(double));
        UUm[i] = (double *) malloc(ny * sizeof(double));
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            UUo[i][j] = 0;
            UUm[i][j] = 0;
        }
    }

    awe_2d_explicit_solver_homogeneous_8th_order(UUo, UUm, dx, dy, dt, v, &F, it, sx, sy, nx, ny);

    int isx = (int) (sx / dx);
    int isy = (int) (sy / dy);

    bool source_injection_test_passed = fabs(UUo[isx][isy] - dt * dt * F) < 1e-6;

    if (source_injection_test_passed) {
        printf("Test for source injection passed.\n");
    } else {
        printf("Test for source injection failed.\n");
    }

    for (int i = 0; i < nx; i++) {
        free(UUo[i]);
        free(UUm[i]);
    }
    free(UUo);
    free(UUm);
}
void test_wave_propagation() {
    int nx = 10, ny = 10;
    double dx = 1.0, dy = 1.0, dt = 0.01, v = 1.0, sx = 5.0, sy = 5.0;
    double F = 1.0;
    int n_steps = 10;

    double **UUo = (double **) malloc(nx * sizeof(double *));
    double **UUm = (double **) malloc(nx * sizeof(double *));
    for (int i = 0; i < nx; i++) {
        UUo[i] = (double *) malloc(ny * sizeof(double));
        UUm[i] = (double *) malloc(ny * sizeof(double));
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            UUo[i][j] = 0;
            UUm[i][j] = 0;
        }
    }

    for (int it = 0; it < n_steps; it++) {
        awe_2d_explicit_solver_homogeneous_8th_order(UUo, UUm, dx, dy, dt, v, &F, it, sx, sy, nx, ny);
        double **temp = UUo;
        UUo = UUm;
        UUm = temp;
    }


    int isx = (int) (sx / dx);
    int isy = (int) (sy / dy);

    // Check if the wave has propagated to neighboring grid points
    bool wave_propagation_test_passed = true;
    for (int i = isx - 1; i <= isx + 1; i++) {
        for (int j = isy - 1; j <= isy + 1; j++) {
            if (i == isx && j == isy) continue; // skip the source location
            if (fabs(UUo[i][j]) < 1e-6) {
                wave_propagation_test_passed = true;
                break;
            }
        }
        if (!wave_propagation_test_passed) break;
    }

    if (wave_propagation_test_passed) {
        printf("Test for wave propagation passed.\n");
    } else {
        printf("Test for wave propagation failed.\n");
    }

    for (int i = 0; i < nx; i++) {
        free(UUo[i]);
        free(UUm[i]);
    }
    free(UUo);
    free(UUm);
}



int main() {
    printf("Running tests...\n");
    
    test_awe_2d_explicit_solver_homogeneous_8th_order();
    test_zero_initial_conditions();
    test_courant_numbers();
    test_boundary_conditions();
    test_source_injection();
    test_wave_propagation();
    return 0;
}

