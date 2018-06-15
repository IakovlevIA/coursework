
#ifndef RENEWAL_OF_COORDINATES_H_
#define RENEWAL_OF_COORDINATES_H_
__device__ vector periodic_boundaries(vector X, double max_valuex, double max_valuey, double max_valuez);

__global__ void New_coordinates(vector *X, vector *P, vector *A, int n,
		int nz, double lattice_constant, double dt,
		double max_valuex, double max_valuey, double max_valuez)
{
        int ix = blockIdx.x * blockDim.x + threadIdx.x;
        int iy = blockIdx.y * blockDim.y + threadIdx.y;
        int iz = blockIdx.z * blockDim.z + threadIdx.z;
        if (ix < n && iy < n && iz < nz){
        	X[iz * n * n + iy * n + ix].x = X[iz * n * n + iy * n + ix].x+
        			P[iz * n * n + iy * n + ix].x*dt+A[iz * n * n + iy * n + ix].x/2*pow(dt,2);
        	X[iz * n * n + iy * n + ix].y = X[iz * n * n + iy * n + ix].y+
        			P[iz * n * n + iy * n + ix].y*dt+A[iz * n * n + iy * n + ix].y/2*pow(dt,2);
        	X[iz * n * n + iy * n + ix].z = X[iz * n * n + iy * n + ix].z+
        			P[iz * n * n + iy * n + ix].z*dt+A[iz * n * n + iy * n + ix].z/2*pow(dt,2);
        	X[iz * n * n + iy * n + ix] = periodic_boundaries(X[iz * n * n + iy * n + ix],
        			max_valuex, max_valuey, max_valuez);
        }
}


#endif /* RENEWAL_OF_COORDINATES_H_ */
