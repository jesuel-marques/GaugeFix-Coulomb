#include <stdio.h>

#include <fields.h>
#include <flags.h>
#include <four_potential.h>
#include <lattice.h>
#include <SU3_ops.h>
#include <types.h>


void SU3_calculate_A(      Mtrx3x3 * restrict U, 
                     const PosVec position, 
                     const LorentzIdx mu,
                           Mtrx3x3 * restrict A) {
    /*  Calculates the vector A_mu(n) field
     The formula is A_mu(n)=((U - U_dagger)/2i)|traceless */
    Mtrx3x3 U_minus_Udagger_trless;

    subtraction_herm_conj_trless_3x3(get_link(U, position, mu), &U_minus_Udagger_trless);
    mult_by_scalar_3x3(-0.5 * I, &U_minus_Udagger_trless, A);
    
}


void SU3_divergence_A(      Mtrx3x3 * restrict U, 
                      const PosVec position, 
                            Mtrx3x3 * restrict div_A) {
    /* Calculates the divergence of the field A on the lattice
    and returns it in div_A. */
    Mtrx3x3 A1;
    Mtrx3x3 A2;

    Mtrx3x3 term_divA;

    set_null_3x3(div_A);
    LorentzIdx mu;
    LOOP_LORENTZ_SPATIAL(mu){
        
        SU3_calculate_A(U,               position,      mu, &A1);
        SU3_calculate_A(U, hop_pos_minus(position, mu), mu, &A2);

        subtraction_3x3(&A1, &A2, &term_divA);

        accumulate_3x3(&term_divA, div_A);  //	Sum of terms over all directions mu.
    }

}