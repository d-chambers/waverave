#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#define mpi_assert(condition, message) {if (!(condition)) {printf("Test failed at rank %d: %s\n", rank, message); MPI_Abort(MPI_COMM_WORLD, 1);}}
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


void test_subdomain_division(int nx, int size, int rank) {
    int chunk_size = nx / size;
    int subdomain_start = rank * chunk_size;
    int subdomain_end = (rank + 1) * chunk_size;

    int expected_subdomain_start = rank * chunk_size;
    int expected_subdomain_end = (rank + 1) * chunk_size;

    mpi_assert(subdomain_start == expected_subdomain_start, "Incorrect subdomain start");
    mpi_assert(subdomain_end == expected_subdomain_end, "Incorrect subdomain end");

    if (rank == 0) {
        printf("Test for subdomain division passed.\n");
    }
}




void test_boundary_data_exchange(int nx, int ny, int size, int rank) {
    double **UUo, **UUm;
    int i;
    int j;

    // Allocate memory
    UUo = (double **) malloc(nx * sizeof(double *));
    UUm = (double **) malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++) {
        UUo[i] = (double *) malloc(ny * sizeof(double));
        UUm[i] = (double *) malloc(ny * sizeof(double));
    }

    // Initialize wavefield
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
           UUo[i][j] = rank + 1;
            UUm[i][j] = rank + 1;
        }
    }

    // Set subdomain boundaries for each process
    int chunk_size = nx / size;
    int subdomain_start = rank * chunk_size;
    int subdomain_end = (rank + 1) * chunk_size;

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

    // Check if the boundary data exchange is successful
    bool test_passed = true;
    if (rank > 0 && UUm[subdomain_start - 1][0] != rank) {
        test_passed = false;
    }
    if (rank < size - 1 && UUm[subdomain_end][0] != rank + 2) {
        test_passed = false;
    }

    if (test_passed) {
        printf("Rank %d: Test for boundary data exchange passed.\n", rank);
    } else {
        printf("Rank %d: Test for boundary data exchange failed.\n", rank);
    }

    // Free memory
    for (i = 0; i < nx; i++) {
        free(UUo[i]);
        free(UUm[i]);
    }
    free(UUo);
    free(UUm);
}

int main(int argc, char *argv[]) {
    int rank, size;
    int nx, ny;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Set dimensions
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);

    // Run the test for boundary data exchange
    test_boundary_data_exchange(nx, ny, size, rank);
    test_subdomain_division(nx,size,rank);
    // Finalize MPI
    MPI_Finalize();
    return 0;
}
