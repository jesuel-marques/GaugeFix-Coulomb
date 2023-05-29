#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <SU3_ops.h>

void calculateCorrelator(Scalar* inverse[4][Nc], double* correlator, const char* type);

#endif  // SPECTRUM_H