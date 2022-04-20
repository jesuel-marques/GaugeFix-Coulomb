#include <complex.h>  //	Standard header C files
#include <stdio.h>
#include <stdlib.h>

#include "../SU3_parameters.h"  //	Simulation parameters
#include "lattice.h"            //	Initialization functions and
                                //	calculations of positions and
                                //	links on the lattice.

#include "SU3_ops.h"  //	SU(3) operations

void SU3_calculate_A(matrix_3x3_double *U, const pos_vec position, const lorentz_index mu, matrix_3x3_double *A) {
    //	Calculates the vector A_mu(n) field
    // The formula is A_mu(n)=((U - U_dagger)/2i)|traceless

    matrix_3x3_double *local_U = get_link(U, position, mu);
    matrix_3x3_double U_dagger;

    SU3_hermitean_conjugate(local_U, &U_dagger);

    matrix_3x3_double U_minus_Udagger;

    //	(U - U_dagger)/2i
    subtraction_3x3(local_U, &U_dagger, &U_minus_Udagger);
    multiplication_by_scalar_3x3(-0.5 * I, &U_minus_Udagger, A);

    //  Subtract the trace part
    matrix_3x3_double trace_part;
    
    set_to_identity_3x3(&trace_part);
    substitution_multiplication_by_scalar_3x3(-( 1.0 / Nc) * trace_3x3(A), &trace_part);

    accumulate_3x3(&trace_part, A);

}

void SU3_divergence_A(matrix_3x3_double *U, const pos_vec position, matrix_3x3_double *div_A) {
    //	Calculates the divergence of the field A on the lattice
    //  and returns it in div_A.

    matrix_3x3_double A1;
    matrix_3x3_double A2;

    matrix_3x3_double term_divA;

    set_to_null_3x3(div_A);
    for (lorentz_index mu = 0; mu < DIM - 1 ; mu++) {
        SU3_calculate_A(U, position, mu, &A1);
        SU3_calculate_A(U, hop_position_negative(position, mu), mu, &A2);

        subtraction_3x3(&A1, &A2, &term_divA);

        accumulate_3x3(&term_divA, div_A);  //	Sum of terms over all directions mu.
    }

}