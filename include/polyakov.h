#ifndef POLYAKOV_H
#define POLYAKOV_H

#include <SU3_ops.h>
#include <geometry.h>

Scalar averagePolyakovLoop(Mtrx3x3* U);

void applyCenterTransformation(Mtrx3x3* U, PosIndex t, CenterElement z);

#endif