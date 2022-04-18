#ifndef SU2OPS_H
#define SU2OPS_H

typedef unsigned short SU2_color_index;

typedef struct {
    double complex m[4];
} matrix_2x2_ck;

void SU2_copy(const matrix_2x2_ck * u, matrix_2x2_ck* u_copy);

void SU2_set_to_null(matrix_2x2_ck * u);
void SU2_set_to_identity(matrix_2x2_ck * u);

void SU2_accumulate(const matrix_2x2_ck * u, matrix_2x2_ck* acc);

void SU2_subtraction(const matrix_2x2_ck* u, const matrix_2x2_ck* v, matrix_2x2_ck* u_minus_v);

inline double SU2_determinant(const matrix_2x2_ck* u);

void SU2_hermitean_conjugate(const matrix_2x2_ck* u, matrix_2x2_ck* u_dagger);

void SU2_multiplication_by_scalar(const matrix_2x2_ck* u, const double alpha, matrix_2x2_ck* alpha_times_u);

void SU2_product(const matrix_2x2_ck* u, const matrix_2x2_ck* v, matrix_2x2_ck* uv);

void SU2_product_three(const matrix_2x2_ck* u, const matrix_2x2_ck* v, const matrix_2x2_ck* w, matrix_2x2_ck* uvw);

void SU2_product_four(const matrix_2x2_ck* u, const matrix_2x2_ck* v, const matrix_2x2_ck* w, const matrix_2x2_ck* x, matrix_2x2_ck* uvwx);

inline void SU2_projection(matrix_2x2_ck* u);

#endif