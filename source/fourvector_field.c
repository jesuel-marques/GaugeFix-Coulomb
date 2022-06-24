#include <complex.h>  //	Standard header C files
#include <stdio.h>
#include <stdlib.h>

#include "../SU3_parameters.h"  //	Simulation parameters
#include "lattice.h"            //	Initialization functions and
                                //	calculations of positions and
                                //	links on the lattice.

#include "SU3_ops.h"  //	SU(3) operations

void SU3_calculate_A(mtrx_3x3 * restrict U, const pos_vec position, const lorentz_idx mu,
                                                                     mtrx_3x3 * restrict A) {
    //	Calculates the vector A_mu(n) field
    // The formula is A_mu(n)=((U - U_dagger)/2i)|traceless

    // mtrx_3x3 *local_U = get_link(U, position, mu);
    // mtrx_3x3 U_dagger;

    // herm_conj_3x3(local_U, &U_dagger);

    mtrx_3x3 U_minus_Udagger_trless;

    //	(U - U_dagger)/2i
    // subtraction_3x3(local_U, &U_dagger, &U_minus_Udagger);
    subtraction_herm_conj_traceless_3x3(get_link(U, position, mu), &U_minus_Udagger_trless);
    mult_by_scalar_3x3(-0.5 * I, &U_minus_Udagger_trless, A);

    //  Subtract the trace part
    // mtrx_3x3 trace_part;
    
    // set_identity_3x3(&trace_part);
    // subst_mult_scalar_3x3(-( 1.0 / Nc) * trace_3x3(A), &trace_part);

    // accumulate_3x3(&trace_part, A);

}

void SU3_divergence_A(mtrx_3x3 * restrict U, const pos_vec position, 
                                            mtrx_3x3 * restrict div_A) {
    //	Calculates the divergence of the field A on the lattice
    //  and returns it in div_A.

    mtrx_3x3 A1;
    mtrx_3x3 A2;

    mtrx_3x3 term_divA;

    set_null_3x3(div_A);
    for (lorentz_idx mu = 0; mu < DIM - 1 ; mu++) {
        SU3_calculate_A(U, position, mu, &A1);
        SU3_calculate_A(U, hop_position_negative(position, mu), mu, &A2);

        subtraction_3x3(&A1, &A2, &term_divA);

        accumulate_3x3(&term_divA, div_A);  //	Sum of terms over all directions mu.
    }

}