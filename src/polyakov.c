#include <SU3_ops.h>
#include <fields.h>
#include <geometry.h>
#include <polyakov.h>
#include <stdio.h>

Scalar polyakovLoop(PosVec position, Mtrx3x3* U) {
    Mtrx3x3 loop;
    setIdentity3x3(&loop);
    LOOP_TEMPORAL(position.pos[T_INDX]) {
        accumRightProd3x3(&loop, getLink(U, position, T_INDX));
    }

    return trace3x3(&loop);
}

Scalar averagePolyakovLoop(Mtrx3x3* U) {
    Scalar average_trace = 0.0;
    Scalar tr_polyakov_loop;

#pragma omp parallel for collapse(DIM - 1) reduction(+ : average_trace)
    for (int k = 0; k < lattice_param.n_SPC; k++) {
        for (int j = 0; j < lattice_param.n_SPC; j++) {
            for (int i = 0; i < lattice_param.n_SPC; i++) {
                tr_polyakov_loop = polyakovLoop(assignPosition(i, j, k, 0), U) / Nc;
                // printf("%lf+I*(%lf)\n", creal(tr_polyakov_loop), cimag(tr_polyakov_loop));
                average_trace += tr_polyakov_loop;
            }
        }
    }

    average_trace /= lattice_param.spatial_volume;
    return average_trace;
}

void applyCenterTransformation(Mtrx3x3* U, PosIndex t, CenterElement z) {
    Mtrx3x3 center_element;

    setIdentity3x3(&center_element);
    substMultScalar3x3(z == ONE                 ? 1.0
                       : z == TWO_PI_OVER_THREE ? cexp(2.0 * I * M_PI / 3.0)
                                                : cexp(4.0 * I * M_PI / 3.0),
                       &center_element);

#pragma omp parallel for collapse(DIM - 1)
    for (int k = 0; k < lattice_param.n_SPC; k++) {
        for (int j = 0; j < lattice_param.n_SPC; j++) {
            for (int i = 0; i < lattice_param.n_SPC; i++) {
                accumLeftProd3x3(&center_element, getLink(U, assignPosition(i, j, k, t), T_INDX));
            }
        }
    }
}