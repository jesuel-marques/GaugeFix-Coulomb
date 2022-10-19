#include <stdio.h>  //	Standard header files in C
#include <stdlib.h>
#include <tgmath.h>

#include <types.h>
#include <math_ops.h>  //	Math operations
#include <SU2_ops.h>
#include <SU3_ops.h>

#include <gauge_fixing.h> // CHANGE THIS. HERE ONLY TO USE SU3_update_sub_Los_Alamos

void print_matrix_3x3(const mtrx_3x3 * restrict u, const char *name, const unsigned short decimal_places) {
    // Prints the matrix on screen with a given number of decimal places and 
    // adds a name on the top

    printf("\n\n %s \n", name);

    printf("{");
    SU3_color_idx a, b;
    
    LOOP_COLOR_3(a){
        printf("{");

        LOOP_COLOR_3(b){
            printf("%.*lf + I*(%.*lf)", decimal_places, creal(u -> m[ELM3x3(a, b)]), 
                                        decimal_places, cimag(u -> m[ELM3x3(a, b)]));
            
            b != Nc -1 ?  printf(", ") : 0 ;
        }

        a != Nc - 1 ?  printf("},\n") : 0 ;
    }

    printf("}}\n\n");

    getchar();
}

void copy_3x3(const mtrx_3x3 * restrict u, mtrx_3x3 * restrict u_copy) {
    // Copies u to u_copy

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){

        u_copy -> m[ELM3x3(a, b)] = u -> m[ELM3x3(a, b)];

    }
}

inline void set_null_3x3(mtrx_3x3 * restrict u) {
    // Sets u to be the 3x3 null matrix

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){

        u -> m[ELM3x3(a, b)] = 0.0;

    }
}

inline void set_identity_3x3(mtrx_3x3 * restrict u) {
    // Sets u to be the identity matrix in SU(3)

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){

                u -> m[ELM3x3(a, b)] = a != b ? 0.0 : 1.0 ;

    }
}

inline void accumulate_3x3(const mtrx_3x3 * restrict u, mtrx_3x3 * restrict acc) {
    // Accumulates the value of u into acc

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){

            acc -> m[ELM3x3(a, b)] += u -> m[ELM3x3(a, b)];
        
    }
}



void subtraction_3x3(const mtrx_3x3 * restrict u, 
                     const mtrx_3x3 * restrict v,
                           mtrx_3x3 * restrict u_minus_v) {
    //  Calculates the difference between matrix u and matrix v
    //  and returns result in u_minus_v

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){

            u_minus_v -> m[ELM3x3(a, b)] =   u -> m[ELM3x3(a, b)] 
                                           - v -> m[ELM3x3(a, b)];

    }
}

scalar trace_3x3(const mtrx_3x3 * restrict u) {
    //	Calculates the trace of the matrix u
    //  and returns result (complex) in tr

    return   u -> m[ELM3x3(0, 0)] 
           + u -> m[ELM3x3(1, 1)] 
           + u -> m[ELM3x3(2, 2)];
}

inline scalar determinant_3x3(const mtrx_3x3 * restrict u) {
    //  Calculates the determinant of the matrix u
    //  and returns result, a complex number, in det.

    //  The determinant is calculated in the usual manner
    //  using Leibniz formula

    return    u -> m[ELM3x3(0, 0)] * (  u -> m[ELM3x3(1, 1)] * u -> m[ELM3x3(2, 2)] 
                                      - u -> m[ELM3x3(1, 2)] * u -> m[ELM3x3(2, 1)]) 
            + u -> m[ELM3x3(0, 1)] * (  u -> m[ELM3x3(1, 2)] * u -> m[ELM3x3(2, 0)] 
                                      - u -> m[ELM3x3(1, 0)] * u -> m[ELM3x3(2, 2)]) 
            + u -> m[ELM3x3(0, 2)] * (  u -> m[ELM3x3(1, 0)] * u -> m[ELM3x3(2, 1)] 
                                      - u -> m[ELM3x3(1, 1)] * u -> m[ELM3x3(2, 0)]);
}

inline void herm_conj_3x3(const mtrx_3x3 * restrict u, 
                                mtrx_3x3 * restrict u_dagger) {
    // Calculates the hermitean conjugate to u
    // and returns result in u_dagger.

    u_dagger -> m[ELM3x3(0, 0)] = conj(u -> m[ELM3x3(0, 0)]);
    u_dagger -> m[ELM3x3(1, 1)] = conj(u -> m[ELM3x3(1, 1)]);
    u_dagger -> m[ELM3x3(2, 2)] = conj(u -> m[ELM3x3(2, 2)]);

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){

            u_dagger -> m[ELM3x3(b, a)] = conj(u -> m[ELM3x3(a, b)]);
            u_dagger -> m[ELM3x3(a, b)] = conj(u -> m[ELM3x3(b, a)]);

    }

}

inline void subtraction_herm_conj_traceless_3x3(const mtrx_3x3 * restrict u, mtrx_3x3 * restrict u_minus_udag_trless) {
    
    scalar tr;
    tr  = u_minus_udag_trless -> m[ELM3x3(0, 0)] = 2 * I * cimag(u -> m[ELM3x3(0, 0)]) ;
    tr += u_minus_udag_trless -> m[ELM3x3(1, 1)] = 2 * I * cimag(u -> m[ELM3x3(1, 1)]) ;
    tr += u_minus_udag_trless -> m[ELM3x3(2, 2)] = 2 * I * cimag(u -> m[ELM3x3(2, 2)]) ;

    tr /= Nc;

    u_minus_udag_trless -> m[ELM3x3(0, 0)] -= tr;
    u_minus_udag_trless -> m[ELM3x3(1, 1)] -= tr;
    u_minus_udag_trless -> m[ELM3x3(2, 2)] -= tr;

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){
            
        u_minus_udag_trless -> m[ELM3x3(a, b)] =       u -> m[ELM3x3(a, b)] 
                                                - conj(u -> m[ELM3x3(b, a)]);
        u_minus_udag_trless -> m[ELM3x3(b, a)] =       u -> m[ELM3x3(b, a)] 
                                                    - conj(u -> m[ELM3x3(a, b)]);
        
    }


}

inline void mult_by_scalar_3x3(const scalar alpha, const mtrx_3x3 * restrict u,
                                            mtrx_3x3 * restrict alpha_times_u) {
    //  Calculates multiplication of 3x3 matrix u by scalar alpha
    //  and returns result in alphatimesu.

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){

            alpha_times_u -> m[ELM3x3(a, b)] = alpha * u -> m[ELM3x3(a, b)];
            //	Mutiplying each entry.

    }
}

inline void subst_mult_scalar_3x3(const scalar alpha, 
                                        mtrx_3x3 * restrict u) {
    //  Calculates multiplication of 3x3 matrix u by scalar alpha
    //  and returns result in u.

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){

        u -> m[ELM3x3(a, b)] *= alpha;
        //	Mutiplying each entry.
        
    }

}

inline void prod_3x3(const mtrx_3x3 * restrict u, 
                     const mtrx_3x3 * restrict v, 
                           mtrx_3x3 * restrict uv) {
    // Calculates product of 2 3x3 matrices u e v
    // and returns result in uv.

    SU3_color_idx a, b, c;
    LOOP_COLOR_3x3(a, b){
        uv -> m[ELM3x3(a, b)] = 0.0;
        LOOP_COLOR_3(c){

            uv -> m[ELM3x3(a, b)] +=   (u -> m[ELM3x3(a, c)]) 
                                     * (v -> m[ELM3x3(c, b)]);
            //  Usual matrix multiplication.
        }
    }
}


void prod_three_3x3(const mtrx_3x3 * restrict u, 
                    const mtrx_3x3 * restrict v,
                    const mtrx_3x3 * restrict w,
                          mtrx_3x3 * restrict uvw) {
    //  Calculates product of 3 3x3 matrices u, v and w
    //  and returns result in uvw.

    SU3_color_idx a, b, c, d;
    LOOP_COLOR_3x3(a, b){
        uvw -> m[ELM3x3(a, b)] = 0.0;
        LOOP_COLOR_3x3(c, d){    

                uvw -> m[ELM3x3(a, b)] +=   (u -> m[ELM3x3(a, c)]) 
                                          * (v -> m[ELM3x3(c, d)]) 
                                          * (w -> m[ELM3x3(d, b)]);
                //  Usual matrix multiplication.
            }
        }

}

void prod_four_3x3(const mtrx_3x3 * restrict u, 
                   const mtrx_3x3 * restrict v,
                   const mtrx_3x3 * restrict w,
                   const mtrx_3x3 * restrict x, 
                         mtrx_3x3 * restrict uvwx){
    //  Calculates product of 3 3x3 matrices u, v and w
    //  and returns result in uvw.

    SU3_color_idx a, b, c, d, e;
    LOOP_COLOR_3x3(a, b){
        uvwx -> m[ELM3x3(a, b)] = 0.0;
        LOOP_COLOR_3(c){
            LOOP_COLOR_3(d){
                LOOP_COLOR_3(e){

                    uvwx -> m[ELM3x3(a, b)] +=   (u -> m[ELM3x3(a, c)]) 
                                                * (v -> m[ELM3x3(c, d)]) 
                                                * (w -> m[ELM3x3(d, e)]) 
                                                * (x -> m[ELM3x3(e, b)]);
                //  Usual matrix multiplication.
                }
            }
        }
    }
}

inline void accum_left_prod_3x3(const mtrx_3x3 * restrict g, mtrx_3x3 * restrict acc_prod) {
    //	Calculates matrix product between g and acc_prod and accumulates result in acc_prod
    mtrx_3x3 aux_prod;
 
    copy_3x3(acc_prod, &aux_prod);
    prod_3x3(g, &aux_prod, acc_prod);

}

inline void accum_right_prod_3x3(mtrx_3x3 * restrict acc_prod, const mtrx_3x3 * restrict g) {
    //	Calculates matrix product between acc_prod and g and accumulates result in acc_prod
    mtrx_3x3 aux_prod;

    copy_3x3(acc_prod, &aux_prod);
    prod_3x3(&aux_prod, g, acc_prod);

}

inline void accum_prod_SU2_3x3(const mtrx_2x2_ck * restrict x_ck, mtrx_3x3 * restrict g, SU3_color_idx a, SU3_color_idx b ){
    //  Accumulates in g the product of a SU(2) Cabbibo-Marinari 
    //  submatrix in Cayley-Klein form x_ck with g
    //  The Cabbibo-Marinari have corresponding line a and column b in 3x3 notation.

    scalar xg1, xg2;
    //  Auxiliary variables
    
    mtrx_2x2 x;
    //  the Cayley-Klein matrix will be first converted to a regular 2x2 matrix
    
    convert_from_ck(x_ck, &x);
    SU3_color_idx c;
    LOOP_COLOR_3(c){

        //  Modifying only the terms in g that will actually receive contributions

        xg1 =   x.m[ELM2x2(0, 0)] * g -> m[ELM3x3(a, c)] 
              + x.m[ELM2x2(0, 1)] * g -> m[ELM3x3(b, c)];
        xg2 =   x.m[ELM2x2(1, 0)] * g -> m[ELM3x3(a, c)] 
              + x.m[ELM2x2(1, 1)] * g -> m[ELM3x3(b, c)];

        g -> m[ELM3x3(a, c)] = xg1;
        g -> m[ELM3x3(b, c)] = xg2;
    }

}

inline void prod_vuwdagger_3x3(const mtrx_3x3 * restrict v, const mtrx_3x3 * restrict u, const mtrx_3x3 * restrict w, mtrx_3x3 * restrict vuwdagger ){
     //	Calculates matrix product between v, u and w dagger.
     
    SU3_color_idx a, b, c, d;
    LOOP_COLOR_3x3(a, b){
        vuwdagger -> m[ELM3x3(a, b)] = 0.0;
        LOOP_COLOR_3x3(c, d){

            vuwdagger -> m[ELM3x3(a, b)] +=       (v -> m[ELM3x3(a, c)]) 
                                            *     (u -> m[ELM3x3(c, d)]) 
                                            * conj(w -> m[ELM3x3(b, d)]);
                                            //  wdagger_eb is conj(w_be)
            
        }
    }
}

inline short power_3x3_binomial(mtrx_3x3 * restrict A, const double omega, mtrx_3x3 * restrict A_to_omega ){
    //  Calculates A^omega according to the formula (which works well for A close to identity)
    // A^omega = Proj_SU3(I + omega * (A-I) + ...)

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){

            A_to_omega -> m[ELM3x3(a, b)] = a == b ? 
                                            1.0 + omega * (A -> m[ELM3x3(a, b)]- 1.0):
                                                  omega * (A -> m[ELM3x3(a, b)]     );
            
    }

    return projection_SU3(A_to_omega);
}


inline double inverse_3x3(const mtrx_3x3 * restrict a, mtrx_3x3 * restrict a_inv){

	//	Calculates the inverse of a 3 by 3 matrix x explicitly
    //  and returns result in a_inv.

	scalar a_det = determinant_3x3(a);
    
    if(a_det != 0.0) {

        a_inv -> m[ELM3x3(0, 0)] =  (a -> m[ELM3x3(1, 1)] * a -> m[ELM3x3(2, 2)] 
                                   - a -> m[ELM3x3(1, 2)] * a -> m[ELM3x3(2, 1)]) / a_det;
        a_inv -> m[ELM3x3(0, 1)] =  (a -> m[ELM3x3(0, 2)] * a -> m[ELM3x3(2, 1)] 
                                   - a -> m[ELM3x3(0, 1)] * a -> m[ELM3x3(2, 2)]) / a_det;
        a_inv -> m[ELM3x3(0, 2)] =  (a -> m[ELM3x3(0, 1)] * a -> m[ELM3x3(1, 2)] 
                                   - a -> m[ELM3x3(0, 2)] * a -> m[ELM3x3(1, 1)]) / a_det;

        a_inv -> m[ELM3x3(1, 0)] =  (a -> m[ELM3x3(1, 2)] * a -> m[ELM3x3(2, 0)] 
                                   - a -> m[ELM3x3(1, 0)] * a -> m[ELM3x3(2, 2)]) / a_det;
        a_inv -> m[ELM3x3(1, 1)] =  (a -> m[ELM3x3(0, 0)] * a -> m[ELM3x3(2, 2)] 
                                   - a -> m[ELM3x3(0, 2)] * a -> m[ELM3x3(2, 0)]) / a_det;
        a_inv -> m[ELM3x3(1, 2)] =  (a -> m[ELM3x3(0, 2)] * a -> m[ELM3x3(1, 0)] 
                                   - a -> m[ELM3x3(0, 0)] * a -> m[ELM3x3(1, 2)]) / a_det;

        a_inv -> m[ELM3x3(2, 0)] =  (a -> m[ELM3x3(1, 0)] * a -> m[ELM3x3(2, 1)] 
                                   - a -> m[ELM3x3(1, 1)] * a -> m[ELM3x3(2, 0)]) / a_det;
        a_inv -> m[ELM3x3(2, 1)] =  (a -> m[ELM3x3(0, 1)] * a -> m[ELM3x3(2, 0)] 
                                   - a -> m[ELM3x3(0, 0)] * a -> m[ELM3x3(2, 1)]) / a_det;
        a_inv -> m[ELM3x3(2, 2)] =  (a -> m[ELM3x3(0, 0)] * a -> m[ELM3x3(1, 1)] 
                                   - a -> m[ELM3x3(0, 1)] * a -> m[ELM3x3(1, 0)]) / a_det;
    
    }
    
    return a_det;

}


inline static int normalize_3_vec(color_3_vec * restrict v){

    //  Takes a color vector and normalizes it
    //  returning result in the same vector

    double r = 0.0;

    SU3_color_idx a;

    LOOP_COLOR_3(a){
        //  Calculates the sum of the absolute
        //  value squared of each component
        r += POW2(fabs(v -> m[a]));

    }
    if ( r != 0 ){
        //  If vector is not null, then normalizes it
        r = 1.0 / sqrt(r);

        LOOP_COLOR_3(a){

            v -> m[a] *= r;

        }
        return EXIT_SUCCESS;
    }
    else{
        return EXIT_FAILURE;
    }

}

inline static void cross_product_conj(color_3_vec u, color_3_vec v, color_3_vec* u_cross_v_star){
    //  Calculates the cross product of u and v and then the 
    //  conjugate of the result and returns result in u_cross_v_star

    //  First component:
    u_cross_v_star -> m[0] = conj(u.m[1] * v.m[2] 
                                - u.m[2] * v.m[1]);

    //  Second component:
    u_cross_v_star -> m[1] = conj(u.m[2] * v.m[0] 
                                - u.m[0] * v.m[2]);

    //  Third component:
    u_cross_v_star -> m[2] = conj(u.m[0] * v.m[1] 
                                - u.m[1] * v.m[0]);

}

inline short projection_SU3(mtrx_3x3 * restrict x) {
    //	Projects matrix x to the group SU(3) 
    //  returning SU(3) matrix x at the end.
    //	Follows method found in openqcd fastsum code.

    short exit_status;

    color_3_vec v0, v1, v2;

    SU3_color_idx b;
    LOOP_COLOR_3(b){

        v0.m[b] = x -> m[ELM3x3(0, b)];
        v1.m[b] = x -> m[ELM3x3(1, b)];

    }

    exit_status = normalize_3_vec(&v0);
    cross_product_conj(v0, v1, &v2);

    exit_status |= normalize_3_vec(&v2);
    cross_product_conj(v2, v0, &v1);

    LOOP_COLOR_3(b){

        x -> m[ELM3x3(0, b)] = v0.m[b];
        x -> m[ELM3x3(1, b)] = v1.m[b];
        x -> m[ELM3x3(2, b)] = v2.m[b];

    }
        
    if(exit_status==EXIT_SUCCESS){
        return EXIT_SUCCESS;
    }
    else{
        return -1;
    }
}

void decompose_algebra_SU3(const mtrx_3x3 * restrict a, mtrx_SU3_alg * restrict a_components) {
    //  Decomposes A in the components of the alfebra of SU(3)
    //  using the Gell-Mann matrices a basis and following the conventions
    //  of Gattringer. The components are then returned in a_components

    // Formulas for the coefficients obtained in Mathematica

    a_components -> m[0] =       0.0;

    a_components -> m[1] =       ( creal(a -> m[ELM3x3(0, 1)]) 
                                 + creal(a -> m[ELM3x3(1, 0)]));

    a_components -> m[2] =       (-cimag(a -> m[ELM3x3(0, 1)]) 
                                 + cimag(a -> m[ELM3x3(1, 0)]));

    a_components -> m[3] =       ( creal(a -> m[ELM3x3(0, 0)]) 
                                 - creal(a -> m[ELM3x3(1, 1)]));

    a_components -> m[4] =       ( creal(a -> m[ELM3x3(0, 2)]) 
                                 + creal(a -> m[ELM3x3(2, 0)]));

    a_components -> m[5] =       (-cimag(a -> m[ELM3x3(0, 2)]) 
                                 + cimag(a -> m[ELM3x3(2, 0)]));

    a_components -> m[6] =       ( creal(a -> m[ELM3x3(1, 2)]) 
                                 + creal(a -> m[ELM3x3(2, 1)]));

    a_components -> m[7] =       (-cimag(a -> m[ELM3x3(1, 2)]) 
                                 + cimag(a -> m[ELM3x3(2, 1)]));

    a_components -> m[8] =       ( creal(a -> m[ELM3x3(0, 0)]) 
                                 + creal(a -> m[ELM3x3(1, 1)]) 
                           - 2.0 * creal(a -> m[ELM3x3(2, 2)])) / sqrt(3);
}

void SU3_CabbiboMarinari_projection(mtrx_3x3 * restrict w, 
                                              mtrx_3x3 * restrict total_update) {
    //	Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
    //	algorithm for SU(3), with a division of the update matrix in submatrices
    //	following the Cabbibo-Marinari trick. Actual update is obtained after a number
    //	of "hits" to be performed one after another.

    mtrx_3x3 w_inv_old;
    mtrx_3x3 total_update_conj;
    //  Calculates the inverse of w in the beginning.
    //  The program will update w successively and to
    //  extract what was the combined update, we can 
    //  multiply from the right by the old inverse.

       
    if(inverse_3x3(w, &w_inv_old)){

        //  Local maximization is attained iteratively in SU(3),
        //  thus we need to make many hits ...
        for (unsigned short hits = 1; hits <= 300; hits++) {

            //	... and each hit contains the Cabbibo-Marinari subdivision
            for (submatrix sub = R; sub <= T; sub++) {
                //	Submatrices are indicated by numbers from 0 to 2
                //  with codenames R, S and T

                SU3_update_sub_LosAlamos(w, sub);
                // projection_SU3(w);
                // printf("re Tr w: %.20lf\n", creal(trace_3x3(w)));
                
            }
        }
    }
    else{
        //  if w has no inverse, update will be given by the identity
        set_identity_3x3(total_update);
        return;
    }
    
    prod_3x3(w, &w_inv_old, &total_update_conj);
    herm_conj_3x3(&total_update_conj, total_update);

    //	Updates matrix to total_update. It is the
    //	accumulated updates from the hits.
}