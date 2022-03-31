#include <complex.h>  //	Standard header C files
#include <stdio.h>
#include <stdlib.h>

#include "../SU3_parameters.h"  //	Simulation parameters
#include "lattice.h"            //	Initialization functions and
                                //	calculations of positions and
                                //	links on the lattice.

#include "SU3_ops.h"  //	SU(3) operations

void SU3_calculate_A(double complex *U, const pos_vec position,const unsigned short mu, double complex *A) {
    //	Calculates the vector A_mu(n) field
    // The formula is A_mu(n)=((U - U_dagger)/2i)|traceless

    double complex *local_U;

    local_U = get_link(U, position, mu);

    double complex *U_dagger = (double complex *)malloc(3 * 3 * sizeof(double complex));
    test_allocation(U_dagger, "SU3_calculate_A");

    SU3_hermitean_conjugate(local_U, U_dagger);

    double complex *U_minus_Udagger = (double complex *)malloc(3 * 3 * sizeof(double complex));
    test_allocation(U_minus_Udagger, "SU3_calculate_A");

    //	(U - U_dagger)/2i
    SU3_subtraction(local_U, U_dagger, U_minus_Udagger);
    SU3_multiplication_by_scalar(-0.5 * I, U_minus_Udagger, A);

    free(U_dagger);
    free(U_minus_Udagger);

    //  Subtract the trace part
    double complex *trace_part = (double complex *)malloc(3 * 3 * sizeof(double complex));
    test_allocation(trace_part, "SU3_calculate_A");
    
    SU3_set_to_identity(trace_part);
    SU3_substitution_multiplication_by_scalar(-(1.0 / 3.0) * SU3_trace(A), trace_part);

    SU3_accumulate(trace_part, A);

    free(trace_part);
}

void SU3_divergence_A(double complex *U,const pos_vec position, double complex *div_A) {
    //	Calculates the divergence of the field A on the lattice
    //  and returns it in div_A.

    double complex *A1 = (double complex *)malloc(3 * 3 * sizeof(double complex));
    test_allocation(A1, "SU3_divergence_A");

    double complex *A2 = (double complex *)malloc(3 * 3 * sizeof(double complex));
    test_allocation(A2, "SU3_divergence_A");

    double complex *term_divA = (double complex *)malloc(3 * 3 * sizeof(double complex));  //  Each term in the sum over all directions   
    test_allocation(term_divA, "SU3_divergence_A");

    SU3_set_to_null(div_A);
    for (unsigned short mu = 1; mu < d; mu++) {
        SU3_calculate_A(U, position, mu, A1);
        SU3_calculate_A(U, hop_position_negative(position, mu), mu, A2);

        SU3_subtraction(A1, A2, term_divA);

        SU3_accumulate(term_divA, div_A);  //	Sum of terms over all directions mu.
    }

    free(A1);
    free(A2);
    free(term_divA);
}