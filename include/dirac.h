#ifndef DIRAC_H

#define DIRAC_H
#include <geometry.h>
#include <stdio.h>
#include <types.h>

#define LOOP_DIRAC(alpha) for (alpha = 0; alpha < 4; alpha++)

typedef struct {
    Scalar m[4][4];
} DiracMatrix;

double invertDiracOperator(double kappa, Mtrx3x3 *U, Scalar *source, Scalar *inverse_column, double tolerance, double (*inversion_algorithm)(void (*)(Scalar *, Scalar *), Scalar *, Scalar *, double, size_t));

void printDiracOP(Mtrx3x3 *U, double kappa, FILE *file_dirac_op);
void initializePointSource(PosVec source_position, int dirac_index, int color_index, double complex *source);

void printInverse(char *inverse_filename, Scalar *inverse[4][3]);

#endif