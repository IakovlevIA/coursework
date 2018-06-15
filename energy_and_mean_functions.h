/*
 * energy_and_mean_functions.h
 *
 *  Created on: 15 июня 2018 г.
 *      Author: deamonde
 */

#ifndef ENERGY_AND_MEAN_FUNCTIONS_H_
#define ENERGY_AND_MEAN_FUNCTIONS_H_

__device__ double min_distance(vector X1, vector X2,
		double max_valuex, double max_valuey, double max_valuez);

__global__ void Energy_and_mean_functions(vector *X, vector *P,
		vector *init_X, vector *init_P, double *kin_E,
		double *pot_E, double *R, double *V, int n,
		int nz, double lattice_constant, double max_valuex,
		double max_valuey, double max_valuez, double epsilon, double sigma)
{
		int ix = blockIdx.x * blockDim.x + threadIdx.x;
        int iy = blockIdx.y * blockDim.y + threadIdx.y;
        int iz = blockIdx.z * blockDim.z + threadIdx.z;
        if (ix < n && iy < n && iz < nz){
        	double r, Z;
        	pot_E[iz * n * n + iy * n + ix] = 0;
        	kin_E[iz * n * n + iy * n + ix] = (pow(P[iz * n * n + iy * n + ix].x,2) +
        			pow(P[iz * n * n + iy * n + ix].y,2) + pow(P[iz * n * n + iy * n + ix].z,2)) / 2;
        	V[iz * n * n + iy * n + ix] = P[iz * n * n + iy * n + ix].x * init_P[iz * n * n + iy * n + ix].x +
						P[iz * n * n + iy * n + ix].y * init_P[iz * n * n + iy * n + ix].y +
						P[iz * n * n + iy * n + ix].z * init_P[iz * n * n + iy * n + ix].z;
        	R[iz * n * n + iy * n + ix] = min_distance(X[iz * n * n + iy * n + ix],
        			init_X[iz * n * n + iy * n + ix], max_valuex, max_valuey, max_valuez);
        	for (int i = 0; i < nz * n * n; i++){
        		if (i != iz * n * n + iy * n + ix){
        			r = min_distance(X[iz * n * n + iy * n + ix], X[i], max_valuex, max_valuey, max_valuez);
        			Z = pow(sigma / r, 6);
        			pot_E[iz * n * n + iy * n + ix] += 4 * epsilon * (pow(Z,2) - Z);
        		}
        	}
        }
}


#endif /* ENERGY_AND_MEAN_FUNCTIONS_H_ */
