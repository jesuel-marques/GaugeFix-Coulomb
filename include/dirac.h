#ifndef DIRAC_H

#define DIRAC_H
#include <types.h>

#define LOOP_DIRAC(alpha) for (alpha = 0; alpha < 4; alpha++)

typedef struct {
    Scalar m[4][4];
} DiracMatrix;

#endif