/*
    four_potential contains routines to calculate the gauge four-potential field and
    related quantities.

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
#include <four_potential.h>
#include <geometry.h>
#include <stdio.h>
#include <types.h>

/* Calculates the vector A_mu(n) field.
   The formula is A_mu=((U - U_dagger)/2i)|traceless */
void calculateA(Mtrx3x3* restrict U,
                const PosVec position,
                const LorentzIdx mu,
                Mtrx3x3* restrict A) {
    /*
     * Calls:
     * =====
     * subtractionHermConjTrless3x3, multByScalar3x3,
     * getLink,
     *
     * Macros:
     * ======
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * U: 	    SU(3) gluon field,
     * PosVec position:     specific position for the calculation,
     * LorentzIdx mu:       specific Lorentz spacetime direction index,
     * Mtrx3x3 * A:         A field calculated.
     *
     * Returns:
     * =======
     *
     */

    Mtrx3x3 U_minus_Udagger_trless;

    subtractionHermConjTrless3x3(getLink(U, position, mu), &U_minus_Udagger_trless);
    multByScalar3x3(-0.5 * I, &U_minus_Udagger_trless, A);
}

/* Calculates the divergence of the field A in a specific position.  */
void DivergenceA(Mtrx3x3* restrict U,
                 const PosVec position,
                 Mtrx3x3* restrict div_A,
                 DivergenceType divergence_type) {
    /*
     * Calls:
     * =====
     * setNull3x3, subtraction3x3, accumulate3x3, substMultScalar3x3
     * getNeighbour,
     * calculateA.
     *
     * Macros:
     * ======
     * LOOP_LORENTZ_SPATIAL
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * U: 	               SU(3) gluon field,
     * PosVec position:                specific position for the calculation,
     * Mtrx3x3 * div_A:                calculated divergence of A at given position.
     * DivergenceType divergence_type: type of divergence to be calculated.
     *
     * Returns:
     * =======
     *
     */

    Mtrx3x3 A1;
    Mtrx3x3 A2;

    Mtrx3x3 term_divA;

    setNull3x3(div_A);
    LorentzIdx mu;

    LOOP_LORENTZ(mu) {
        if (mu != T_INDX || divergence_type == QUADRI) {
            calculateA(U, position, mu, &A1);
            calculateA(U, getNeighbour(position, mu, REAR), mu, &A2);
            subtraction3x3(&A1, &A2, &term_divA);
            if (mu == T_INDX) {
                substMultScalar3x3(lattice_param.divergence_anisotropy, &term_divA);
            }
            accumulate3x3(&term_divA, div_A);
        }
    }
}