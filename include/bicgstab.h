#ifndef BICGSTAB_H
#define BICGSTAB_H

void initializePointSource(PosVec source_position, int dirac_index, int color_index, double complex *source);
void BiCGStab(int dirac_index, int color_index, double kappa, double tolerance, Mtrx3x3 *U, Scalar *source, Scalar *inverse_column);

#endif  // BICGSTAB_H