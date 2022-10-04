#include <tgmath.h>
#include <stdio.h>  //	Standard header files in C
#include <stdlib.h>

#include <math_ops.h>  //	Math operations
#include <lattice.h>
#include <SU2_ops.h>
#include <SU3_ops.h>

void print_matrix_3x3(const mtrx_3x3 * restrict u, const char *name, const unsigned short decimal_places) {
    // Prints the matrix on screen with a given number of decimal places and 
    // adds a name on the top

    printf("\n\n %s \n", name);

    printf("{");
    for (SU3_color_idx  a = 0; a < Nc; a++) {
        printf("{");

        for (SU3_color_idx  b = 0; b < Nc; b++) {
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

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            u_copy -> m[ELM3x3(a, b)] = u -> m[ELM3x3(a, b)];

        }
    }
}

void convert_in_work_cfg_3x3(const in_cfg_data_type * restrict u_in, 
                                work_mtrx_data_type * restrict u_work) {

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            u_work -> m[ELM3x3(a, b)] = (work_data_type) u_in -> m[ELM3x3(a, b)];

        }
    }
}

void convert_work_out_cfg_3x3(const work_mtrx_data_type * restrict u_work, 
                                      out_cfg_data_type * restrict u_out) {


    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            u_out -> m[ELM3x3(a, b)] = (out_data_type) u_work -> m[ELM3x3(a, b)];

        }
    }
}


void convert_in_work_gt_3x3(const in_gt_data_type * restrict g_in, 
                              work_mtrx_data_type * restrict g_work) {

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            g_work -> m[ELM3x3(a, b)] = (work_data_type) g_in -> m[ELM3x3(a, b)];

        }
    }
}

void convert_work_out_gt_3x3(const work_mtrx_data_type * restrict g_work, 
                                      out_gt_data_type * restrict g_out) {


    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            g_out -> m[ELM3x3(a, b)] = (out_data_type) g_work -> m[ELM3x3(a, b)];

        }
    }
}

inline void set_null_3x3(mtrx_3x3 * restrict u) {
    // Sets u to be the 3x3 null matrix

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            u -> m[ELM3x3(a, b)] = 0.0;

        }
    }
}

inline void set_identity_3x3(mtrx_3x3 * restrict u) {
    // Sets u to be the identity matrix in SU(3)

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

                u -> m[ELM3x3(a, b)] = a != b ? 0.0 : 1.0 ;

        }
    }
}

inline void accumulate_3x3(const mtrx_3x3 * restrict u, mtrx_3x3 * restrict acc) {
    // Accumulates the value of u into acc

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            acc -> m[ELM3x3(a, b)] += u -> m[ELM3x3(a, b)];
        
        }
    }
}



void subtraction_3x3(const mtrx_3x3 * restrict u, 
                     const mtrx_3x3 * restrict v,
                           mtrx_3x3 * restrict u_minus_v) {
    //  Calculates the difference between matrix u and matrix v
    //  and returns result in u_minus_v

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            u_minus_v -> m[ELM3x3(a, b)] =   u -> m[ELM3x3(a, b)] 
                                           - v -> m[ELM3x3(a, b)];

        }
    }
}

work_data_type trace_3x3(const mtrx_3x3 * restrict u) {
    //	Calculates the trace of the matrix u
    //  and returns result (complex) in tr

    return   u -> m[ELM3x3(0, 0)] 
           + u -> m[ELM3x3(1, 1)] 
           + u -> m[ELM3x3(2, 2)];
}

inline work_data_type determinant_3x3(const mtrx_3x3 * restrict u) {
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

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < a; b++) {

            u_dagger -> m[ELM3x3(b, a)] = conj(u -> m[ELM3x3(a, b)]);
            u_dagger -> m[ELM3x3(a, b)] = conj(u -> m[ELM3x3(b, a)]);

        }
    }
}

inline void subtraction_herm_conj_traceless_3x3(const mtrx_3x3 * restrict u, mtrx_3x3 * restrict u_minus_udag_trless) {
    
    work_data_type tr;
    tr  = u_minus_udag_trless -> m[ELM3x3(0, 0)] = 2 * I * cimag(u -> m[ELM3x3(0, 0)]) ;
    tr += u_minus_udag_trless -> m[ELM3x3(1, 1)] = 2 * I * cimag(u -> m[ELM3x3(1, 1)]) ;
    tr += u_minus_udag_trless -> m[ELM3x3(2, 2)] = 2 * I * cimag(u -> m[ELM3x3(2, 2)]) ;

    tr /= Nc;

    u_minus_udag_trless -> m[ELM3x3(0, 0)] -= tr;
    u_minus_udag_trless -> m[ELM3x3(1, 1)] -= tr;
    u_minus_udag_trless -> m[ELM3x3(2, 2)] -= tr;

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < a; b++) {
            
            u_minus_udag_trless -> m[ELM3x3(a, b)] =       u -> m[ELM3x3(a, b)] 
                                                    - conj(u -> m[ELM3x3(b, a)]);
            u_minus_udag_trless -> m[ELM3x3(b, a)] =       u -> m[ELM3x3(b, a)] 
                                                    - conj(u -> m[ELM3x3(a, b)]);
        
        }
    }


}

inline void mult_by_scalar_3x3(const work_data_type alpha, const mtrx_3x3 * restrict u,
                                            mtrx_3x3 * restrict alpha_times_u) {
    //  Calculates multiplication of 3x3 matrix u by scalar alpha
    //  and returns result in alphatimesu.

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            alpha_times_u -> m[ELM3x3(a, b)] = alpha * u -> m[ELM3x3(a, b)];
            //	Mutiplying each entry.

        }
    }
}

inline void subst_mult_scalar_3x3(const work_data_type alpha, 
                                                                mtrx_3x3 * restrict u) {
    //  Calculates multiplication of 3x3 matrix u by scalar alpha
    //  and returns result in u.

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            u -> m[ELM3x3(a, b)] *= alpha;
            //	Mutiplying each entry.
            
        }
    }
}

inline void prod_3x3(const mtrx_3x3 * restrict u, 
                     const mtrx_3x3 * restrict v, 
                           mtrx_3x3 * restrict uv) {
    // Calculates product of 2 3x3 matrices u e v
    // and returns result in uv.

    for (SU3_color_idx  a = 0; a < Nc; a++) {          //  lines
        for (SU3_color_idx  b = 0; b < Nc; b++) {      //  columns
            uv -> m[ELM3x3(a, b)] = 0.0;
            for (SU3_color_idx  c = 0; c < Nc; c++) {  //  dummy index

                uv -> m[ELM3x3(a, b)] +=   (u -> m[ELM3x3(a, c)]) 
                                         * (v -> m[ELM3x3(c, b)]);
                //  Usual matrix multiplication.
            }
        }
    }
}

void prod_three_3x3(const mtrx_3x3 * restrict u, 
                    const mtrx_3x3 * restrict v,
                    const mtrx_3x3 * restrict w,
                          mtrx_3x3 * restrict uvw) {
    //  Calculates product of 3 3x3 matrices u, v and w
    //  and returns result in uvw.

    for (SU3_color_idx  a = 0; a < Nc; a++) {              //  lines
        for (SU3_color_idx  b = 0; b < Nc; b++) {          //  columns
            uvw -> m[ELM3x3(a, b)] = 0.0;
            for (SU3_color_idx  c = 0; c < Nc; c++) {      //  dummy index 1
                for (SU3_color_idx  e = 0; e < Nc; e++) {  //  dummy index 2

                    uvw -> m[ELM3x3(a, b)] +=   (u -> m[ELM3x3(a, c)]) 
                                              * (v -> m[ELM3x3(c, e)]) 
                                              * (w -> m[ELM3x3(e, b)]);
                    //  Usual matrix multiplication.
                }
            }
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

    for (SU3_color_idx  a = 0; a < Nc; a++) {              //  lines
        for (SU3_color_idx  b = 0; b < Nc; b++) {          //  columns
            uvwx -> m[ELM3x3(a, b)] = 0.0;
            for (SU3_color_idx  c = 0; c < Nc; c++) {            //  dummy index 1
                for (SU3_color_idx  e = 0; e < Nc; e++) {        //  dummy index 2
                    for (SU3_color_idx  f = 0; f < Nc; f++) {    //  dummy index 3

                        uvwx -> m[ELM3x3(a, b)] +=   (u -> m[ELM3x3(a, c)]) 
                                                   * (v -> m[ELM3x3(c, e)]) 
                                                   * (w -> m[ELM3x3(e, f)]) 
                                                   * (x -> m[ELM3x3(f, b)]);
                    //  Usual matrix multiplication.
                    }
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

    work_data_type xg1, xg2;
    //  Auxiliary variables
    
    mtrx_2x2 x;
    //  the Cayley-Klein matrix will be first converted to a regular 2x2 matrix
    
    convert_from_ck(x_ck, &x);
   
    for(SU3_color_idx c = 0 ; c < Nc ; c++){

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
     
     for (SU3_color_idx  a = 0; a < Nc; a++) {              //  lines
        for (SU3_color_idx  b = 0; b < Nc; b++) {          //  columns
            vuwdagger -> m[ELM3x3(a, b)] = 0.0;
            for (SU3_color_idx  c = 0; c < Nc; c++) {      //  dummy index 1
                for (SU3_color_idx  e = 0; e < Nc; e++) {  //  dummy index 2

                    vuwdagger -> m[ELM3x3(a, b)] +=       (v -> m[ELM3x3(a, c)]) 
                                                    *     (u -> m[ELM3x3(c, e)]) 
                                                    * conj(w -> m[ELM3x3(b, e)]);
                                                    //  wdagger_eb is conj(w_be)
                }
            }
        }
    }
}

inline short power_3x3_binomial(mtrx_3x3 * restrict A, const double omega, mtrx_3x3 * restrict A_to_omega ){
    //  Calculates A^omega according to the formula (which works well for A close to identity)
    // A^omega = Proj_SU3(I + omega * (A-I) + ...)

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            A_to_omega -> m[ELM3x3(a, b)] = a == b ? 
                                            1.0 + omega * (A -> m[ELM3x3(a, b)]- 1.0):
                                                  omega * (A -> m[ELM3x3(a, b)]     );
            
        }
    }

    return projection_SU3(A_to_omega);
}


inline double inverse_3x3(const mtrx_3x3 * restrict a, mtrx_3x3 * restrict a_inv){

	//	Calculates the inverse of a 3 by 3 matrix x explicitly
    //  and returns result in a_inv.

	work_data_type a_det = determinant_3x3(a);
    
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

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        //  Calculates the sum of the absolute
        //  value squared of each component
        r += POW2(fabs(v -> m[a]));

    }
    if ( r != 0 ){
        //  If vector is not null, then normalizes it
        r = 1.0 / sqrt(r);

        for (SU3_color_idx  a = 0; a < Nc; a++) {

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

    for (SU3_color_idx  b = 0; b < Nc; b++) {

        v0.m[b] = x -> m[ELM3x3(0, b)];
        v1.m[b] = x -> m[ELM3x3(1, b)];

    }

    exit_status = normalize_3_vec(&v0);
    cross_product_conj(v0, v1, &v2);

    exit_status |= normalize_3_vec(&v2);
    cross_product_conj(v2, v0, &v1);

    for (SU3_color_idx  b = 0; b < Nc; b++) {

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