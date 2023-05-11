/*
    fields takes care of field operations on a lattice.

    Copyright (C) 2023  Jesuel Marques

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Contact: jesuel.leal@usp.br

 */

#include <../settings.h>
#include <SU3_ops.h>
#include <fields.h>
#include <geometry.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <types.h>

/* Allocated the memory for 3x3 with a given number of elements. */
Mtrx3x3 *allocate3x3Field(unsigned elements) {
    /*
     * Calls:
     * =====
     * calloc.
     *
     * Macros:
     * ======
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * unsigned elements:   number of 3x3 matrices of the field to be allocated.
     *
     * Returns:
     * =======
     * A pointer to the allocated field. This pointer can be NULL, if the internal
     * calloc finds some error which prevents it from allocating the necessary memory.
     *
     */

    const Mtrx3x3 *field = (Mtrx3x3 *)calloc(elements, sizeof(Mtrx3x3));

    return field;
}

/* Sets all 3x3 matrices entries in field to unit matrices. */
int setFieldToIdentity(const Mtrx3x3 *restrict field, unsigned elements) {
    /*
     * Calls:
     * =====
     * setIdentity3x3.
     *
     * Macros:
     * ======
     * OMP_PARALLEL_FOR.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * field:     field composed as a collection of 3x3 matrices,
     * unsigned elements:   number of 3x3 matrices of the field.
     *
     * Returns:
     * =======
     * 0 if successful.
     * Integer error-code otherwise:
     *      -2: the field passed was a NULL pointer.
     *
     */

    if (field == NULL) {
        return -2;
    }

#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
    for (int i = 0; i < elements; i++)
        setIdentity3x3(field + i);

    return 0;
}

/* Sets all 3x3 matrices entries in field to special unitary random matrices using de
ranlux random number generator. */
int setFieldSU3Random(const Mtrx3x3 *restrict field, unsigned elements) {
    /*
     * Calls:
     * =====
     * setSU3Random.
     *
     * Macros:
     * ======
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * field:     field composed as a collection of 3x3 matrices,
     * unsigned elements:   number of 3x3 matrices of the field.
     *
     * Returns:
     * =======
     * 0 if successful.
     * Integer error-code otherwise:
     *      -2: the field passed was a NULL pointer.
     *
     */

    if (field == NULL) {
        return -2;
    }

    for (int i = 0; i < elements; i++)
        setSU3Random(field + i);

    return 0;
}

/* Copies a field to field_copy. */
int copyField(const Mtrx3x3 *restrict field,
              unsigned elements,
              const Mtrx3x3 *restrict field_copy) {
    /*
     * Calls:
     * =====
     * copy3x3.
     *
     * Macros:
     * ======
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * field:         field composed as a collection of 3x3 matrices,
     * unsigned elements:       number of 3x3 matrices of the field,
     * Mtrx3x3 * field_copy:    copy of field.
     *
     * Returns:
     * =======
     * 0 if successful.
     * Integer error-code otherwise:
     *      -2: the field or field_copy passed was a NULL pointer.
     *
     */

    if (field == NULL || field_copy == NULL) {
        return -2;
    }

#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
    for (int i = 0; i < elements; i++)
        copy3x3(field + i, field_copy + i);

    return 0;
}

/* Reunitarizes all the matrices in the field by projecting them all to SU(3). */
int reunitarizeField(const Mtrx3x3 *restrict field, unsigned elements) {
    /*
     * Calls:
     * =====
     * projectSU3.
     *
     * Macros:
     * ======
     * NUM_THREADS.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * field:     field composed as a collection of 3x3 matrices,
     * unsigned elements:   number of 3x3 matrices of the field.
     *
     * Returns:
     * =======
     * 0 if successful.
     * Integer error-code otherwise:
     *      -1: the projection didn't work, probably because there is a
     *          null matrix somewhere.
     *      -2: the field passed was a NULL pointer.
     *
     */

    if (field == NULL) {
        return -2;
    }

    int exit_status = 0;
#pragma omp parallel for num_threads(NUM_THREADS) \
    reduction(| : exit_status) schedule(dynamic)
    for (int i = 0; i < elements; i++)
        exit_status |= projectSU3(field + i);

    if (exit_status != 0) {
        return -1;
    }
    return exit_status;
}

/* Gives the average determinant for a field of 3x3 matrices. */
Scalar averageFieldDet(const Mtrx3x3 *restrict field, unsigned elements) {
    /*
     * Calls:
     * =====
     * determinant3x3.
     *
     * Macros:
     * ======
     * NUM_THREADS.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * field:     field composed as a collection of 3x3 matrices,
     * unsigned elements:   number of 3x3 matrices of the field.
     *
     * Returns:
     * =======
     * The lattice averaged determinant for the field if successful.
     * Integer error-code otherwise:
     *      -2: the field passed was a NULL pointer.
     *
     */

    if (field == NULL) {
        return -2;
    }
    Scalar det = 0.0;

#pragma omp parallel for num_threads(NUM_THREADS) \
    schedule(dynamic) reduction(+ : det)
    for (int i = 0; i < elements; i++)
        det += determinant3x3(field + i);

    det /= (Scalar)elements;

    return det;
}

/* Gets a pointer to a link at the given position and mu. */
Mtrx3x3 *getLink(const Mtrx3x3 *restrict U,
                 const PosVec position,
                 const LorentzIdx mu) {
    /*
     * Calls:
     * =====
     * if CHECK_POSITION_BOUNDS set
     * validPositionQ.
     *
     * Macros:
     * ======
     * CHECK_POSITION_BOUNDS, GET_LINK.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * U:         SU(3) gluon field,
     * PosVec position:     position at which the link variable is requested,
     * LorentzIdx mu:       Lorentz direction for which the link variable is requested.
     *
     * Returns:
     * =======
     * The gauge link at that particular point and mu as a pointer to a 3x3 matrix.
     *
     */

#ifdef CHECK_POSITION_BOUNDS
    if (positionmuValidQ(position, mu))
#endif  // CHECK_POSITION_BOUNDS
        return GET_LINK(U, position, mu);
#ifdef CHECK_POSITION_BOUNDS
    else {
        return NULL;
    }
#endif  // CHECK_POSITION_BOUNDS
}

/* Gets forward or backward link at given position and mu. */
void getLinkMatrix(const Mtrx3x3 *restrict U,
                   const PosVec position,
                   const LorentzIdx mu,
                   Direction dir,
                   const Mtrx3x3 *restrict u) {
    /*
     * Calls:
     * =====
     * getNeighbour,
     * copy3x3, hermConj3x3,
     * getLink.
     *
     * Macros:
     * ======
     * CHECK_POSITION_BOUNDS.
     *
     * if CHECK_POSITION_BOUNDS set
     * positionmuValidQ.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * U:	        SU(3) gluon field,
     * PosVec position:     position at which the link is required,
     * LorentzIdx mu:       Lorentz direction at which the link is required,
     * Direction dir:       FRONT if the forward link REAR if the backward link,
     * Mtrx3x3 * u:         the link with coordinates position and mu.
     *
     * Returns:
     * =======
     *
     */

#ifdef CHECK_POSITION_BOUNDS
    if (positionmuValidQ(position, mu)) {
#endif  // CHECK_POSITION_BOUNDS
        if (dir == FRONT) {
            copy3x3(getLink(U, position, mu), u);
            //	Link in the positive way is what is stored in U

        } else if (dir == REAR) {
            hermConj3x3(getLink(U, getNeighbour(position, mu, REAR), mu), u);
            //	U_(-mu)(n) = (U_mu(n - mu))^\dagger
        } else {
            u = NULL;
        }
#ifdef CHECK_POSITION_BOUNDS
    } else {
        u = NULL;
    }
#endif  // CHECK_POSITION_BOUNDS
}

/* Gets a pointer to the gauge transformation at the given position. */
Mtrx3x3 *getGaugetransf(const Mtrx3x3 *restrict G,
                        const PosVec position) {
    /*
     * Calls:
     * =====
     * if CHECK_POSITION_BOUNDS set
     * validPositionQ.
     *
     * Macros:
     * ======
     * CHECK_POSITION_BOUNDS, GET_GT.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * G:         SU(3) gauge transformation field,
     * PosVec position:     position at which the gauge transformation is requested.
     *
     * Returns:
     * =======
     * The gauge transformation at that particular point as a pointer to a 3x3 matrix.
     *
     */

#ifdef CHECK_POSITION_BOUNDS
    if (validPositionQ(position)) {
#endif  // CHECK_POSITION_BOUNDS
        return GET_GT(G, position);
#ifdef CHECK_POSITION_BOUNDS
    } else {
        return NULL;
    }
#endif  // CHECK_POSITION_BOUNDS
}

void applyGaugeTransformationU(Mtrx3x3 *restrict U,
                               Mtrx3x3 *restrict G) {
    /*
     * Description:
     * ===========
     * Applies gauge transformation G to gluon-field U.
     *
     * Calls:
     * =====
     * prod_vuwdagger3x3, copy3x3,
     * getNeighbour,
     * getLink, getGaugetransf.
     *
     * Macros:
     * ======
     * LOOP_TEMPORAL, T_INDX, LOOP_SPATIAL, LOOP_LORENTZ.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * U:	    SU(3) gluon field,
     * Mtrx3x3 * G:	    SU(3) gauge-transformation field.
     *
     * Returns:
     * =======
     *
     */

    PosIndex t;

#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
    LOOP_TEMPORAL(t) {
        Mtrx3x3 *g;
        Mtrx3x3 *u;
        Mtrx3x3 u_updated;

        PosVec position;

        position.pos[T_INDX] = t;

        LOOP_SPATIAL(position) {
            g = getGaugetransf(G, position);
            LorentzIdx mu;
            LOOP_LORENTZ(mu) {
                //	U'_mu(x) = g(x) . U_mu(x) . gdagger(x + mu)
                prod_vuwdagger3x3(g,
                                  u = getLink(U, position, mu),
                                  getGaugetransf(G, getNeighbour(position, mu, FRONT)),
                                  &u_updated);

                copy3x3(&u_updated, u);
            }
        }
    }
}

void applyRandomGaugeTransformationU(Mtrx3x3 *restrict U) {
    /*
     * Description:
     * ===========
     * Applies random gauge transformation to gluon-field U.
     *
     * Calls:
     * =====
     * prod_vuwdagger3x3, copy3x3,
     * getNeighbour,
     * getLink, getGaugetransf.
     *
     * Macros:
     * ======
     * LOOP_TEMPORAL, T_INDX, LOOP_SPATIAL, LOOP_LORENTZ.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * U:	    SU(3) gluon field,
     * Mtrx3x3 * G:	    SU(3) gauge-transformation field.
     *
     * Returns:
     * =======
     *
     */

    PosIndex t;

    Mtrx3x3 *G = allocate3x3Field(lattice_param.volume);
    setFieldSU3Random(G, lattice_param.volume);

    applyGaugeTransformationU(U, G);
    free(G);
}