
#ifndef INITIALIZATION_OF_IMPULSES_H_
#define INITIALIZATION_OF_IMPULSES_H_

__host__ void Init_impulse(vector *P, vector *init_P,
		int n, int nz, double p0)
{
	double last_px = 0, last_py = 0, last_pz = 0;
	double px_rand, py_rand, pz_rand;
	for (int i = 0; i < n * n * nz; i++){
		px_rand = rand();
		py_rand = rand();
		pz_rand = rand();
		P[i].x = 2 * (px_rand / RAND_MAX - 0.5) * p0;
		P[i].y = 2 * (py_rand / RAND_MAX - 0.5) * p0;
		P[i].z = 2 * (pz_rand / RAND_MAX - 0.5) * p0;
		last_px += P[i].x;
		last_py += P[i].y;
		last_pz += P[i].z;
	}

	for (int i = 0; i < n * n * nz; i++){
		P[i].x -= last_px / (n * n * nz);
		P[i].y -= last_py / (n * n * nz);
		P[i].z -= last_pz / (n * n * nz);
		init_P[i].x = P[i].x;
		init_P[i].y = P[i].y;
		init_P[i].z = P[i].z;
	}
}


#endif /* INITIALIZATION_OF_IMPULSES_H_ */
