#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

void awe_2d_explicit_solver_homogeneous_8th_order(double **UUo, double **UUm, double dx, double dy, double dt, double v, double *F, int it, double sx, double sy, int nx, int ny, int subdomain_start, int subdomain_end) {
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
    if (isx >= subdomain_start && isx < subdomain_end) {
        UUo[isx][isy] += dt * dt * F[it];
    }

    // Update solution
    for (i = MAX(4, subdomain_start); i < MIN(nx - 4, subdomain_end); i++) {
        for (j = 4; j < ny - 4; j++) {
            UUm[i][j] = 2 * UUo[i][j] - UUm[i][j] + Cx2 * (
                -1.0 / 560 * UUo[i - 4][j] + 8.0 / 315 * UUo[i - 3][j] - 1.0 / 5 * UUo[i - 2][j] + 8.0 / 5 * UUo[i - 1][j] - 205.0 / 72 * UUo[i][j] + 8.0 / 5 * UUo[i + 1][j] - 1.0 / 5 * UUo[i + 2][j] + 8.0 / 315 * UUo[i + 3][j] - 1.0 / 560 * UUo[i + 4][j]
            ) + Cy2 * (
                -1.0 / 560 * UUo[i][j - 4] + 8.0 / 315 * UUo[i][j - 3] - 1.0 / 5 * UUo[i][j - 2] + 8.0 / 5 * UUo[i][j - 1] - 205.0 / 72 * UUo[i][j] + 8.0 / 5 * UUo[i][j + 1] - 1.0 / 5 * UUo[i][j + 2] + 8.0 / 315 * UUo[i][j + 3] - 1.0 / 560 * UUo[i][j + 4]
            );
        }
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    double cc;
    int Lx, Ly, nx, ny, it, nt;
    double **UUo, **UUm;
    double *F;
    double dx, dy, dt, v, sx, sy;
    int subdomain_start, subdomain_end;
    double start_time, end_time;
    int i,j;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Set parameters
    Lx = 500;
    Ly = 300;
    nx = atoi (argv[1]);
    ny = atoi (argv[2]);
    nt = 600;
    dx = Lx / nx;
    dy = Ly / ny;
    cc = 0.5;
    dt = cc * dx /v;
    v = 1500.0;
    sx = nx * dx / 2.0;
    sy = ny * dy / 2.0;

    // Allocate memory
    UUo = (double **) malloc(nx * sizeof(double *));
    UUm = (double **) malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++) {
        UUo[i] = (double *) malloc(ny * sizeof(double));
        UUm[i] = (double *) malloc(ny * sizeof(double));
    }
    F = (double *) malloc(nt * sizeof(double));

    // Initialize wavefield and source function
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            UUo[i][j] = 0.0;
            UUm[i][j] = 0.0;
        }
    }
    for (i = 0; i < nt; i++) {
        F[i] = exp(-pow((i - nt / 2) * dt, 2) / (pow(0.1, 2)));
    }

    // Set subdomain boundaries for each process
    int chunk_size = nx / size;
    subdomain_start = rank * chunk_size;
    subdomain_end = (rank + 1) * chunk_size;

    // Start timer
    start_time = MPI_Wtime();

    // Run simulation
    for (it = 0; it < nt; it++) {
        awe_2d_explicit_solver_homogeneous_8th_order(UUo, UUm, dx, dy, dt, v, F, it, sx, sy, nx, ny, subdomain_start, subdomain_end);

        // Exchange data between subdomains
        MPI_Request req;
        if (rank > 0) {
            MPI_Isend(UUm[subdomain_start], ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req);
            MPI_Recv(UUm[subdomain_start - 1], ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Isend(UUm[subdomain_end - 1], ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &req);
            MPI_Recv(UUm[subdomain_end], ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Swap UUm and UUo pointers
        double **temp = UUm;
        UUm = UUo;
        UUo = temp;
    }

    // End timer
    end_time = MPI_Wtime();

    // Print execution time for each process
    if (rank == 0) {
        printf("Total execution time (seconds): %f\n", end_time - start_time);
    }

    // Clean up and finalize
    for (i = 0; i < nx; i++) {
        free(UUo[i]);
        free(UUm[i]);
    }
    free(UUo);
    free(UUm);
    free(F);

    MPI_Finalize();
    return 0;
}


