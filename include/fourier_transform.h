#ifndef FOURIER_TRANSFORM_H
#define FOURIER_TRANSFORM_H

#include <geometry.h>

DiracColorMatrix fourier_transform(DiracColorMatrix* restrict x_space_vector, double momentum[DIM]);

#endif
