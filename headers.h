
#ifndef HEADERS_H_
#define HEADERS_H_

#define BLOCK_SIZE 8
#define BLOCK_SIZE_Z 2

typedef struct vec {
    double x;
    double y;
    double z;
    double buf;
} vector;

#include "initialization_of_coordinates.h"
#include "initialization_of_impulses.h"
#include "renewal_of_coordinates.h"
#include "initialization_of_accelerations.h"
#include "renewal_of_impulses_and_accelerations.h"
#include "energy_and_mean_functions.h"
#include "min_distance.h"
#include "new_acceleration.h"
#include "periodic_boundaries.h"

#endif /* HEADERS_H_ */
