/*
 * new_acceleration.h
 *
 *  Created on: 15 июня 2018 г.
 *      Author: deamonde
 */

#ifndef NEW_ACCELERATION_H_
#define NEW_ACCELERATION_H_


__device__ vector new_acceleration(vector X1, vector X2,
		double max_valuex, double max_valuey,
		double max_valuez, double sigma, double epsilon)
{
	double r, rx, ry, rz, Z;
	vector new_a;
	rx = X1.x - X2.x;
	if (abs(rx) > max_valuex / 2) {
		rx = rx / abs(rx) * (abs(rx) - max_valuex);
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
	Z = pow(sigma/r,6);
	new_a.x = 24 * epsilon * Z * rx / r / r * (2 * Z - 1);
	new_a.y = 24 * epsilon * Z * ry / r / r * (2 * Z - 1);
	new_a.z = 24 * epsilon * Z * rz / r / r * (2 * Z - 1);
	return(new_a);
}


#endif /* NEW_ACCELERATION_H_ */
