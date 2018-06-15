/*
 * min_distance.h
 *
 *  Created on: 15 июня 2018 г.
 *      Author: deamonde
 */

#ifndef MIN_DISTANCE_H_
#define MIN_DISTANCE_H_

__device__ double min_distance(vector X1, vector X2,
		double max_valuex, double max_valuey, double max_valuez)
{
	double r, rx, ry, rz;
	rx = X1.x - X2.x;
	if (abs(rx) > max_valuex / 2) {
		rx = rx / abs(rx)*(abs(rx) - max_valuex);
	}
	ry = X1.y - X2.y;
	if (abs(ry) > max_valuey / 2) {
		ry = ry / abs(ry) * (abs(ry) - max_valuey);
	}
	rz = X1.z - X2.z;
	if (abs(rz) > max_valuez / 2) {
		rz = rz / abs(rz) * (abs(rz) - max_valuez);
	}
	r = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));
	return(r);
}



#endif /* MIN_DISTANCE_H_ */
