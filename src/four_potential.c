#include <stdio.h>

#include <fields.h>
#include <flags.h>
#include <four_potential.h>
#include <lattice.h>
#include <SU3_ops.h>
#include <types.h>


void calculateA(Mtrx3x3 * restrict U, 
                const PosVec position, 
                const LorentzIdx mu,
                Mtrx3x3 * restrict A) {
    /*  Calculates the vector A_mu(n) field
     The formula is A_mu(n)=((U - U_dagger)/2i)|traceless */
    Mtrx3x3 U_minus_Udagger_trless;

    subtractionHermConjTrless3x3(getLink(U, position, mu), &U_minus_Udagger_trless);
    multByScalar3x3(-0.5 * I, &U_minus_Udagger_trless, A);
}


void divergenceA(Mtrx3x3 * restrict U, 
                 const PosVec position, 
                 Mtrx3x3 * restrict div_A) {
    /* Calculates the divergence of the field A on the lattice
    and returns it in div_A. */
    Mtrx3x3 A1;
    Mtrx3x3 A2;

    Mtrx3x3 term_divA;

    setNull3x3(div_A);
    LorentzIdx mu;
    LOOP_LORENTZ_SPATIAL(mu) {
        
        calculateA(U,              position,            mu, &A1);
        calculateA(U, getNeighbour(position, mu, REAR), mu, &A2);

        subtraction3x3(&A1, &A2, &term_divA);

        accumulate3x3(&term_divA, div_A);  //	Sum of terms over all directions mu.

    }
}