
#ifndef INITIALIZATION_OF_COORDINATES_H_
#define INITIALIZATION_OF_COORDINATES_H_

__global__ void Init_coordinates(vector *X, vector *init_X, int n,int nz,double lattice_constant)
{
        int ix = blockIdx.x * blockDim.x + threadIdx.x;
        int iy = blockIdx.y * blockDim.y + threadIdx.y;
        int iz = blockIdx.z * blockDim.z + threadIdx.z;
        if (ix < n && iy < n && iz < nz){
        	X[iz * n * n + iy * n + ix].x = lattice_constant * ix;
        	X[iz * n * n + iy * n + ix].y = lattice_constant * iy;
        	X[iz * n * n + iy * n + ix].z = lattice_constant * iz;
        	init_X[iz * n * n + iy * n + ix] = X[iz * n * n + iy * n + ix];
        }
}


#endif /* INITIALIZATION_OF_COORDINATES_H_ */
