#include <dirac.h>
#include <fourier_transform.h>
#include <geometry.h>
#include <stdlib.h>
#include <tgmath.h>

#define DIGITS_MATRIX 5

void setnullDiracColorMatrix(DiracColorMatrix* restrict u) {
    unsigned int a, b;

    LOOP_DCXDC(a, b) {
        u->m[ELEM_12x12(a, b)] = 0.0;
    }
}

void accumulateDiracColorMatrix(DiracColorMatrix* restrict u, DiracColorMatrix* restrict acc) {
    unsigned int a, b;
    LOOP_DCXDC(a, b) {
        acc->m[ELEM_12x12(a, b)] += u->m[ELEM_12x12(a, b)];
    }
}

inline void multByScalarDiracColorMatrix(const Scalar scalar,
                                         DiracColorMatrix* restrict u,
                                         DiracColorMatrix* restrict scalar_times_u) {
    unsigned int a_alpha, b_beta;
    LOOP_DCXDC(a_alpha, b_beta) {
        scalar_times_u->m[ELEM_12x12(a_alpha, b_beta)] = scalar * u->m[ELEM_12x12(a_alpha, b_beta)];
    }
}

inline double ScalarProductPosVecMom(PosVec position, double momentum[DIM]) {
    double scalar_product = 0.0;
    LorentzIdx mu;
    LOOP_LORENTZ(mu) {
        scalar_product += position.pos[mu] * momentum[mu];
    }

    return scalar_product;
}

void printMatrix12x12(DiracColorMatrix* restrict u) {
    unsigned short a, b;
    printf("{");
    LOOP_DC(a) {
        printf("{");
        LOOP_DC(b) {
            printf("%.*lf + I*(%.*lf)", DIGITS_MATRIX, creal(u->m[ELEM_12x12(a, b)]),
                   DIGITS_MATRIX, cimag(u->m[ELEM_12x12(a, b)]));

            b != 4 * Nc - 1 ? printf(", ") : 0;
        }
        a != 4 * Nc - 1 ? printf("},\n") : 0;
    }
    printf("}}\n\n");
}

DiracColorMatrix fourierTransform(DiracColorMatrix* restrict x_space_vector, double momentum[DIM]) {
    PosVec position;

    DiracColorMatrix prop_momentum, term;
    setnullDiracColorMatrix(&prop_momentum);
    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            multByScalarDiracColorMatrix(cexp(-I * ScalarProductPosVecMom(position, momentum)), (ELEM_VEC_POS(x_space_vector, position)), &term);
            accumulateDiracColorMatrix(&term, &prop_momentum);
        }
    }
    // printf("%lf %lf %lf %lf\n", momentum[T_INDX], momentum[X_INDX], momentum[Y_INDX], momentum[Z_INDX]);
    // printMatrix12x12(&prop_momentum);
    // getchar();
    return prop_momentum;
}