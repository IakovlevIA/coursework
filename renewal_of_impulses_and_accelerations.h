
#ifndef RENEWAL_OF_IMPULSES_AND_ACCELERATIONS_H_
#define RENEWAL_OF_IMPULSES_AND_ACCELERATIONS_H_

__device__ vector new_acceleration(vector X1, vector X2,
		double max_valuex, double max_valuey,
		double max_valuez, double sigma, double epsilon);

__global__ void New_impulse_and_acceleration(vector *X, vector *P,
		vector *A, int n, int nz, double lattice_constant, double max_valuex,
		double max_valuey, double max_valuez, double epsilon, double dt, double sigma)
{
        int ix = blockIdx.x * blockDim.x + threadIdx.x;
        int iy = blockIdx.y * blockDim.y + threadIdx.y;
        int iz = blockIdx.z * blockDim.z + threadIdx.z;
        if (ix < n && iy < n && iz < nz){
		double ax_new = 0, ay_new = 0, az_new = 0;
		vector a_new;
		for (int i = 0; i < nz * n * n; i++){
    			if (i != iz*n*n+iy*n+ix){
    				a_new=new_acceleration(X[iz * n * n + iy * n + ix], X[i],
    						max_valuex, max_valuey, max_valuez, sigma, epsilon);
    				ax_new+=a_new.x;
    				ay_new+=a_new.y;
    				az_new+=a_new.z;
			}
		}
		P[iz * n * n + iy * n + ix].x = P[iz * n * n + iy * n + ix].x +
				(ax_new + A[iz * n * n + iy * n + ix].x) / 2 * dt;
		P[iz * n * n + iy * n + ix].y = P[iz * n * n + iy * n + ix].y +
				(ay_new + A[iz * n * n + iy * n + ix].y) / 2 * dt;
		P[iz * n * n + iy * n + ix].z = P[iz * n * n + iy * n + ix].z +
				(az_new + A[iz * n * n + iy * n + ix].z) / 2 * dt;
		A[iz * n * n + iy * n + ix].x=ax_new;
		A[iz * n * n + iy * n + ix].y=ay_new;
		A[iz * n * n + iy * n + ix].z=az_new;
        }
}



#endif /* RENEWAL_OF_IMPULSES_AND_ACCELERATIONS_H_ */
