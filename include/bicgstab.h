#ifndef BICGSTAB_H
#define BICGSTAB_H

#include <geometry.h>

void initializePointSource(PosVec source_position, int dirac_index, int color_index, double complex *source);
void BiCGStab(double kappa, Mtrx3x3 *U, Scalar *source, Scalar *inverse_column, double tolerance);

#endif  // BICGSTAB_H