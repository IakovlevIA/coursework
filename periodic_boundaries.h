
#ifndef PERIODIC_BOUNDARIES_H_
#define PERIODIC_BOUNDARIES_H_

__device__ vector periodic_boundaries(vector X,
		double max_valuex, double max_valuey, double max_valuez)
{
	if (X.x > max_valuex){
		X.x = X.x - max_valuex;
	}
	if (X.x < 0){
		X.x = X.x + max_valuex;
	}
	if (X.y > max_valuey){
		X.y = X.y - max_valuey;
	}
	if (X.y < 0){
		X.y = X.y + max_valuey;
	}
	if (X.z > max_valuez){
		X.z = X.z - max_valuez;
	}
	if (X.z < 0){
		X.z = X.z + max_valuez;
	}
	return(X);
}

#endif /* PERIODIC_BOUNDARIES_H_ */
