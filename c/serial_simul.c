#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

int main() {
    // Define spatial grid
    double Lx = 500, Ly = 300;
    int nx = 801, ny = 481;
    double dx = Lx / (nx - 1), dy = Ly / (ny - 1);
    double x[nx], y[ny];
    int i, j, it;

    for (i = 0; i < nx; i++) {
        x[i] = i * dx;
    }

    for (j = 0; j < ny; j++) {
        y[j] = j * dy;
    }

    double v = 3000; // Velocity (m/s)

    // Allocate memory for wavefields
    double **UUo = (double **) malloc(nx * sizeof(double *));
    double **UUm = (double **) malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++) {
        UUo[i] = (double *) malloc(ny * sizeof(double));
        UUm[i] = (double *) malloc(ny * sizeof(double));
    }

    // Initialize wavefields
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            UUo[i][j] = 0;
            UUm[i][j] = 0;
        }
    }

    // Time stepping parameters
    double CC = 0.5; // Courant #
    int nt = 600;
    double dt = CC * dx / v;

    double t[nt];

    for (it = 0; it < nt; it++) {
        t[it] = it * dt;
    }

    // Define initial waveform
    double peak_freq = 60; // sigma for Ricker wavelet
    double F[nt];

    for (it = 0; it < nt; it++) {
        F[it] = (1 -2* pow(M_PI,2)*pow(peak_freq,2)*pow(t[it],2)) * exp(-pow(M_PI,2)* pow(peak_freq,2)*pow(t[it],2));
    }

    double sx = Lx / 2, sy = Ly / 2;

    // Iterate over solution
    // Save the wavefield to a binary file
    FILE *file = fopen("wavefield_serial.bin", "wb");
    if (file == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }

    for (it = 0; it < nt; it++) {
        awe_2d_explicit_solver_homogeneous_8th_order(UUo, UUm, dx, dy, dt, v, F, it, sx, sy, nx, ny);
        
        for (i = 0; i < nx; i++) {
            fwrite(UUo[i], sizeof(double), ny, file);
        }

        // Swap wavefields for the next iteration
        double **temp = UUm;
        UUm = UUo;
        UUo = temp;
    }

    fclose(file);


    // Free memory
    for (i = 0; i < nx; i++) {
        free(UUo[i]);
        free(UUm[i]);
    }
    free(UUo);
    free(UUm);

    return 0;
}

