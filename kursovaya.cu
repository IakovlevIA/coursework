#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <curand_kernel.h>
#include <time.h>
#include <fstream>
#include "headers.h"

int main(void) {
	FILE *out1, *out2, *inp, *begin, *end;
	inp = fopen("inp.dat", "r");
	out1 = fopen("particles.dat", "w");
	out2 = fopen("param.dat", "w");
	begin = fopen("begin.dat", "w");
	end = fopen("end.dat", "w");
	double dt, p0, lattice_constant,
	epsilon, sigma, T = 0, kb=1,
	kin_E, pot_E, total_E, V, R, Lx,
	max_valuex, max_valuey, max_valuez;
	int n, nz, iter;
	int t0 = time(0);
	vector *h_ParticleX, *d_ParticleX;
	vector *h_ParticleP, *d_ParticleP;
	vector *h_ParticleA, *d_ParticleA;
	vector *d_init_particleX;
	vector *h_init_particleP, *d_init_particleP;
	double *h_R, *d_R;
	double *h_V, *d_V;
	double *h_kin_E, *d_kin_E;
	double *h_pot_E, *d_pot_E;

	//input
	fscanf(inp, "%d", &n);
	fscanf(inp, "%d", &nz);
	fscanf(inp, "%d", &iter);
	fscanf(inp, "%lf", &dt);
	fscanf(inp, "%lf", &Lx);
	fscanf(inp, "%lf", &p0);
	fscanf(inp, "%lf", &epsilon);
	fscanf(inp, "%lf", &sigma);
	fclose(inp);

	size_t vsize = nz * n * n * sizeof(vector);
	size_t size = nz * n * n * sizeof(double);

	// host memory allocation
	h_ParticleX = (vector *) malloc(vsize);
	h_ParticleP = (vector *) malloc(vsize);
	h_ParticleA = (vector *) malloc(vsize);
	h_init_particleP = (vector *) malloc(vsize);
	h_R = (double *) malloc(size);
	h_V = (double *) malloc(size);
	h_kin_E = (double *) malloc(size);
	h_pot_E = (double *) malloc(size);

	// device memory allocation
    cudaMalloc(&d_ParticleX, vsize);
	cudaMalloc(&d_ParticleP, vsize);
	cudaMalloc(&d_ParticleA, vsize);
	cudaMalloc(&d_init_particleX, vsize);
	cudaMalloc(&d_init_particleP, vsize);
	cudaMalloc(&d_R, size);
	cudaMalloc(&d_V, size);
	cudaMalloc(&d_kin_E, size);
	cudaMalloc(&d_pot_E, size);

	srand(time(NULL));

	//initialization of threads geometry
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE_Z);
    double GridSize = (n * n) / BLOCK_SIZE;
    double GridSizeZ = nz / BLOCK_SIZE_Z;
    dim3 dimGrid(GridSize, GridSize, GridSizeZ);

	//system size initialization
	lattice_constant = Lx / (n - 1);
	max_valuex = lattice_constant * n;
	max_valuey = max_valuex;
	max_valuez = lattice_constant * nz;

	Init_coordinates<<< dimGrid, dimBlock>>>(d_ParticleX,
			d_init_particleX, n, nz, lattice_constant);

	Init_impulse(h_ParticleP, h_init_particleP, n, nz, p0);

	cudaMemcpy(d_ParticleP, h_ParticleP, vsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_init_particleP, h_init_particleP, vsize, cudaMemcpyHostToDevice);

	Init_acceleration<<< dimGrid, dimBlock>>>(d_ParticleX, d_ParticleA, n, nz,
			lattice_constant, max_valuex, max_valuey, max_valuez, epsilon, sigma);

	//main part
	for (int time_count = 0; time_count < iter; time_count++){

	//calculation of energy and mean functions
		if (time_count % 10 == 0){
			Energy_and_mean_functions<<< dimGrid, dimBlock>>>(d_ParticleX,
					d_ParticleP, d_init_particleX, d_init_particleP, d_kin_E,
					d_pot_E, d_R, d_V, n, nz, lattice_constant, max_valuex,
					max_valuey, max_valuez, epsilon, sigma);

			cudaMemcpy(h_ParticleX, d_ParticleX, vsize, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_ParticleP, d_ParticleP, vsize, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_ParticleA, d_ParticleA, vsize, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_kin_E, d_kin_E, size, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_pot_E, d_pot_E, size, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_R, d_R, size, cudaMemcpyDeviceToHost);
			cudaMemcpy(h_V, d_V, size, cudaMemcpyDeviceToHost);

			pot_E = 0;
			kin_E = 0;
			V = 0;
			R = 0;
			for(int i = 0; i < nz * n * n; i++){
				kin_E += h_kin_E[i];
				pot_E += h_pot_E[i];
				R += h_R[i];
				V += h_V[i];
			}
			total_E = kin_E + pot_E * 0.5;
			T += kin_E / nz / n / n / kb;
			V = V / (n * n * nz);
			R = R / (n * n * nz);

		//output
			for(int i = 0; i < nz * n * n; i++){
				fprintf(out1,"%10.10lf   %10.10lf  %10.10lf  %10.10lf "
						"  %10.10lf   %10.10lf  %10.10lf  %10.10lf  %10.10lf\n",
					h_ParticleX[i].x, h_ParticleX[i].y, h_ParticleX[i].z,
					h_ParticleP[i].x, h_ParticleP[i].y, h_ParticleP[i].z,
					h_ParticleA[i].x, h_ParticleA[i].y, h_ParticleA[i].z);
			}
			fprintf(out1,"\n \n");
			fprintf(out2,"%d %10.10lf %10.10lf %10.10lf %10.10lf %10.10lf %10.10lf\n",
					time_count, kin_E, pot_E*0.5, total_E, T/time_count, V, R);

			if (time_count == 0){
				for(int i = 0; i < n * n * nz; i++){
					fprintf(begin,"%10.10lf   %10.10lf  %10.10lf\n",
							h_ParticleX[i].x, h_ParticleX[i].y, h_ParticleX[i].z);
				}
			}
			if (time_count == iter - 10){
				for(int i = 0; i < n * n * nz; i++){
					fprintf(end,"%10.10lf   %10.10lf  %10.10lf\n",
							h_ParticleX[i].x, h_ParticleX[i].y, h_ParticleX[i].z);
				}
			}
		}

		New_coordinates<<< dimGrid, dimBlock>>>(d_ParticleX, d_ParticleP, d_ParticleA,
				n, nz, lattice_constant, dt, max_valuex, max_valuey, max_valuez);

		New_impulse_and_acceleration<<< dimGrid, dimBlock>>>(d_ParticleX, d_ParticleP,
				d_ParticleA, n, nz, lattice_constant, max_valuex, max_valuey, max_valuez,
				epsilon, dt, sigma);
	}
	printf("time %d s\n ", time(0) - t0);
	fclose(out1);
	fclose(out2);
	fclose(begin);
	fclose(end);
	cudaFree(d_ParticleA);
	cudaFree(d_ParticleX);
	cudaFree(d_ParticleP);
	cudaFree(d_kin_E);
	cudaFree(d_pot_E);
	cudaFree(d_R);
	cudaFree(d_V);
}
