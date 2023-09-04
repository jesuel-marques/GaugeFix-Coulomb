#ifndef DIRAC_H

#define DIRAC_H
#include <SU3_ops.h>
#include <geometry.h>
#include <stdio.h>
#include <types.h>

typedef unsigned short DiracIdx;

#define LOOP_DIRAC(alpha) for (alpha = 0; alpha < 4; alpha++)

#define ELEM_VEC_POS(v, position) (v) + ((((position.pos[T_INDX]) * lattice_param.n_SPC + (position.pos[Z_INDX])) * lattice_param.n_SPC + (position.pos[Y_INDX])) * lattice_param.n_SPC + (position.pos[X_INDX]))

#define ELEM_VEC_POSDC(v, position, dirac, color) ((v) + (((((position.pos[T_INDX]) * lattice_param.n_SPC + (position.pos[Z_INDX])) * lattice_param.n_SPC + (position.pos[Y_INDX])) * lattice_param.n_SPC + (position.pos[X_INDX])) * 4 + (dirac)) * Nc + (color))

typedef struct DiracMatrix {
    Scalar m[4 * 4];
} DiracMatrix;

typedef struct DiracColorMatrix {
    Scalar m[4 * Nc * 4 * Nc];
} DiracColorMatrix;

typedef struct DiracIdxPair {
    DiracIdx alpha;
    DiracIdx beta;
} DiracIdxPair;

typedef struct SparseArrayEntry {
    size_t index[2];
    Scalar value;
} SparseArrayEntry;

typedef struct SparseArray {
    SparseArrayEntry *entries;
    size_t number_of_entries;
} SparseArray;

#define ELEM_DCxDC(alpha, a, beta, b) ((((alpha) * Nc + (a)) * 4 + (beta)) * Nc + (b))
#define ELEM_12x12(a, b) (a) * 12 + (b)

#define LOOP_DCXDC(a, b)         \
    for (a = 0; a < 4 * Nc; a++) \
        for (b = 0; b < 4 * Nc; b++)

#define LOOP_DC(a) \
    for (a = 0; a < 4 * Nc; a++)

void setupDiracOperator(double kappa, double c_SW, Mtrx3x3 *U);

double invertDiracOperator(Scalar *source, Scalar *inverse_column, const double tolerance, double (*inversion_algorithm)(Scalar *(*)(Scalar *, Scalar *), Scalar *, Scalar *, double, size_t));

void destroyDiracOperator();

void printDiracOperator(FILE *file_dirac_op);

DiracMatrix prodDirac(DiracMatrix m1, DiracMatrix m2);
void initializePointSource(PosVec source_position, DiracIdx dirac_index, MtrxIdx3 color_index, Scalar *source);

void printInverse(const char *restrict inverse_filename, Scalar *restrict inverse[4][Nc]);

DiracMatrix calculateMatrixSlash(double vector[DIM]);

DiracMatrix colortraceFromDiracColor(DiracColorMatrix m_DC);
Scalar diractraceFromDirac(DiracMatrix m_D);

void aK_from_ap(double ap[4], double aK[4]);
void aQ_from_ap(double ap[4], double aQ[4]);

double dotprod(double ap[4], double aq[4]);

void printDiracMatrix(DiracMatrix u);

#endif