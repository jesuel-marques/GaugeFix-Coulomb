#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include <SU3_ops.h>

double averageSpatialPlaquette (Mtrx3x3 * U),
       averageTemporalPlaquette(Mtrx3x3 * U);

#endif  //MEASUREMENT_H