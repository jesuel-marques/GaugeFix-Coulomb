#ifndef FOURIER_TRANSFORM_H
#define FOURIER_TRANSFORM_H

#include <geometry.h>

DiracColorMatrix fourierTransform(DiracColorMatrix* restrict x_space_vector, double momentum[DIM]);

void printMatrix12x12(DiracColorMatrix* restrict u);

#endif
