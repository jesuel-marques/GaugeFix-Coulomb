#ifndef FOURIER_TRANSFORM_H
#define FOURIER_TRANSFORM_H

#include <geometry.h>

DiracColorMatrix fourierTransform(DiracColorMatrix* restrict x_space_vector, double momentum[DIM]);

#endif
