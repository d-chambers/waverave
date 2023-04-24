#include <math.h>
#include <stdio.h>

void awe_2d_explicit_solver_heterogeneous_8th_order(double **UUo, double **UUm, int nx, int ny, double dx, double dy, double dt, double **v, double *F, int it, double sx, double sy) {

//inputs: 
//UUo: Acoustic pressure vector (nx,ny) at time step n
//UUm: Acoustic pressure vector (nx,ny) at time step n-1
//dx : Spatial sampling in x
//dy : Spatial sampling in y        
//dt : Temporal sampling
//v  : Heterogeneous propagation velocity (nx,ny)
//F  : source time function (nt)
//it : Time index
//sx : Location of source in x (meters)
//sy : Location of source in y (meters)


//Calculate Courant numbers (squared)
	double dtdx2 = pow(dt / dx, 2);
	double dtdy2 = pow(dt / dy, 2);

//Source injection
	int isx = (int)(sx/dx);
	int isy = (int)(sy/dy);

//inject wavelet
	UUo[isx][isy] += dt * dt * F[it];

//update solution
	for (int i = 4; i < nx - 4; i++) {
		for (int j = 4; j < ny - 4; j++) {
			UUm[i][j] = 2 * UUo[i][j] - UUm[i][j]
			+ dtdx2 * pow(v[i][j], 2) * (
			-1 / 560.0 * UUo[i - 4][j] +
			8 / 315.0 * UUo[i - 3][j] -
			1 / 5.0 * UUo[i - 2][j] +
			8 / 5.0 * UUo[i - 1][j] -
			205 / 72.0 * UUo[i][j] +
			8 / 5.0 * UUo[i + 1][j] -
			1 / 5.0 * UUo[i + 2][j] +
			8 / 315.0 * UUo[i + 3][j] -
			1 / 560.0 * UUo[i + 4][j]
			)
			+ dtdy2 * pow(v[i][j], 2) * (
			-1 / 560.0 * UUo[i][j - 4] +
			8 / 315.0 * UUo[i][j - 3] -
			1 / 5.0 * UUo[i][j - 2] +
			8 / 5.0 * UUo[i][j - 1] -
			205 / 72.0 * UUo[i][j] +
			8 / 5.0 * UUo[i][j + 1] -
			1 / 5.0 * UUo[i][j + 2] +
			8 / 315.0 * UUo[i][j + 3] -
			1 / 560.0 * UUo[i][j + 4]
			);
		}
	}

//swap pointers
		double **temp;
		temp = UUo;
		UUo = UUm;
		UUm = temp;
}

double **allocate_2d_array(int nx, int ny) {
    double **arr = (double **)malloc(nx * sizeof(double *));
    for (int i = 0; i < nx; i++) {
        arr[i] = (double *)malloc(ny * sizeof(double));
    }
    return arr;
}

void deallocate_2d_array(double **arr, int nx) {
    for (int i = 0; i < nx; i++) {
        free(arr[i]);
    }
    free(arr);
}
