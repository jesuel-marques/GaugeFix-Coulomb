#ifndef SU2OPS_H
#define SU2OPS_H

typedef unsigned short SU2_color_index;

typedef struct {
    double complex m[4];
} matrix_2x2_ck;

void copy_2x2(const matrix_2x2_ck * u, matrix_2x2_ck* u_copy);

void set_to_null_2x2(matrix_2x2_ck * u);
void set_to_identity_2x2(matrix_2x2_ck * u);

void accumulate_2x2(const matrix_2x2_ck * u, matrix_2x2_ck* acc);

void subtraction_2x2(const matrix_2x2_ck* u, const matrix_2x2_ck* v, matrix_2x2_ck* u_minus_v);

inline double determinant_2x2(const matrix_2x2_ck* u);

void hermitean_conjugate_2x2(const matrix_2x2_ck* u, matrix_2x2_ck* u_dagger);

void multiplication_by_scalar_2x2(const matrix_2x2_ck* u, const double alpha, matrix_2x2_ck* alpha_times_u);

void product_2x2(const matrix_2x2_ck* u, const matrix_2x2_ck* v, matrix_2x2_ck* uv);

void product_three_2x2(const matrix_2x2_ck* u, const matrix_2x2_ck* v, const matrix_2x2_ck* w, matrix_2x2_ck* uvw);

void product_four_2x2(const matrix_2x2_ck* u, const matrix_2x2_ck* v, const matrix_2x2_ck* w, const matrix_2x2_ck* x, matrix_2x2_ck* uvwx);

inline void SU2_projection(matrix_2x2_ck* u);

#endif