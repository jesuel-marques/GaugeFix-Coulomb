/*
    measurement performs measurements on lattice quantities like plaquettes.

    Copyright (C) 2023  Jesuel Marques

    This program is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this
    program.  If not, see <https://www.gnu.org/licenses/>.

    Contact: jesuel.leal@usp.br

 */

#include <SU3_ops.h>
#include <fields.h>
#include <geometry.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <types.h>

extern GeometricParameters lattice_param;

/* Calculates the trace of a specified plaquette. */
Scalar TrPlaquette(Mtrx3x3* restrict U,
                   const PosVec position,
                   const LorentzIdx mu,
                   const LorentzIdx nu) {
    /*
     * Calls:
     * =====
     * getNeighbour, getLinkMatrix,
     * prodFour3x3, trace3x3.
     *
     * Macros:
     * ======
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 *  U:        SU(3) gluon field,
     * PosVec position:     position to which the plaquette is attached,
     * LorentzIdx mu:       first direction of the plaquette specification,
     * LorentzIdx nu:       second direction of the plaquette specification.
     *
     * Returns:
     * =======
     * The trace of the plaquette as a double precision complex number.
     *
     */

    Mtrx3x3 plaquette;

    Mtrx3x3 ua, ub, uc, ud;

    const PosVec position_plus_mu = getNeighbour(position, mu, FRONT);

    getLinkMatrix(U, position, mu, FRONT, &ua);
    getLinkMatrix(U, position_plus_mu, nu, FRONT, &ub);
    getLinkMatrix(U, getNeighbour(position_plus_mu, nu, FRONT), mu, REAR, &uc);
    getLinkMatrix(U, getNeighbour(position, nu, FRONT), nu, REAR, &ud);

    prodFour3x3(&ua, &ub, &uc, &ud, &plaquette);

    return trace3x3(&plaquette) / Nc;
}

/* Average plaquette value for any given mu and nu */
Scalar averagePlaquettegeneric(Mtrx3x3* U, LorentzIdx mu, LorentzIdx nu) {
    /*
     * Calls:
     * =====
     * TrPlaquette.
     *
     * Macros:
     * ======
     * LOOP_SPATIAL, T_INDX, NUM_THREADS.
     *
     * Global Variables:
     * ================
     * lattice_param.
     *
     * Parameters:
     * ==========
     * Mtrx3x3 *  U:        SU(3) gluon field,
     * LorentzIdx mu:       first direction of the plaquette specification,
     * LorentzIdx nu:       second direction of the plaquette specification.
     *
     * Returns:
     * =======
     * The average plaquette value for given mu and nu as a double precision complex
     * number.
     *
     */

    Scalar plaq_ave = 0.0;

// Parallelizing by slicing the time extent
#pragma omp parallel for reduction(+ : plaq_ave) \
    num_threads(NUM_THREADS) schedule(dynamic)
    for (PosIndex t = 0; t < lattice_param.n_T; t++) {
        PosVec position;
        position.pos[T_INDX] = t;

        LOOP_SPATIAL(position) {
            plaq_ave += TrPlaquette(U, position, mu, nu);
        }
    }

    plaq_ave /= ((double)lattice_param.volume);

    return plaq_ave;
}

/* Calculates the average plaquette on the lattice of specific type (total,
   spatial, temporal ) */
Scalar averagePlaquette(Mtrx3x3* restrict U, char* type) {
    /*
     * Calls:
     * =====
     * fprintf, strcmp,
     * validGeometricParametersQ,
     * averagePlaquettegeneric.
     *
     * Macros:
     * ======
     * LOOP_LORENTZ_SPATIAL, T_INDX, X_INDX, Y_INDX, Z_INDX.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 *  U:        SU(3) gluon field,
     * char * type:         type of plaquette average requested. Can be either total,
     *                      spatial, temporal. Other values are ignored.
     *
     * Returns:
     * =======
     * The average trace of the plaquettes of the gluon field as a double precision
     * complex number for the specified type.
     *
     */

    if (!validGeometricParametersQ()) {
        fprintf(stderr, "Error in geometric parameters\n");
    }

    double average = 0.0;
    unsigned short count = 0;
    if (!strcmp(type, "total") || !strcmp(type, "temporal")) {
        LorentzIdx i;
        LOOP_LORENTZ_SPATIAL(i) {
            average += averagePlaquettegeneric(U, T_INDX, i);
            count++;
        }
    }
    if (!strcmp(type, "total") || !strcmp(type, "spatial")) {
        LorentzIdx i, j;
        LOOP_LORENTZ_SPATIAL(i) {
            LOOP_LORENTZ_SPATIAL(j) {
                if (i < j) {
                    average += averagePlaquettegeneric(U, i, j);
                    count++;
                }
            }
        }
    }

    return average / (double)count;
}