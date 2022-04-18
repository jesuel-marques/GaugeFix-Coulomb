#ifndef SU3OPS_H
#define SU3OPS_H

typedef unsigned short SU3_color_index;
typedef unsigned short SU3_color_alg_index;

typedef struct {
    double complex m[Nc * Nc];
} matrix_3x3_double;


typedef struct {
    float complex m[Nc * Nc];
} matrix_3x3_float;

void SU3_print_matrix(const matrix_3x3_double *u, const char *name);

void SU3_copy(const matrix_3x3_double *u, matrix_3x3_double *u_copy);

void SU3_convert_fd(const matrix_3x3_float *u_float, matrix_3x3_double *u_double);

void SU3_convert_df(const matrix_3x3_double *u_double, matrix_3x3_float *u_float);

void SU3_set_to_null(matrix_3x3_double *u);

void SU3_set_to_identity(matrix_3x3_double *u);

void SU3_accumulate(const matrix_3x3_double *u, matrix_3x3_double *acc);

void SU3_subtraction(const matrix_3x3_double *u, const matrix_3x3_double *v, matrix_3x3_double *u_minus_v);

inline extern double complex SU3_trace(const matrix_3x3_double *u);

inline extern double complex SU3_determinant(const matrix_3x3_double *u);

inline extern void SU3_hermitean_conjugate(const matrix_3x3_double *u, matrix_3x3_double *u_dagger);

inline extern void SU3_multiplication_by_scalar(const double complex alpha, const matrix_3x3_double *u, matrix_3x3_double *alpha_times_u);
inline extern void SU3_substitution_multiplication_by_scalar(const double complex alpha, matrix_3x3_double *u);

inline extern void SU3_product(const matrix_3x3_double *u, const matrix_3x3_double *v, matrix_3x3_double *uv);

void SU3_product_three(const matrix_3x3_double *u, const matrix_3x3_double *v, const matrix_3x3_double *w, matrix_3x3_double *uvw);

void SU3_product_four(const matrix_3x3_double *u, const matrix_3x3_double *v, const matrix_3x3_double *w, const matrix_3x3_double *x, matrix_3x3_double *uvwx);

inline extern void SU3_accumulate_left_product(const matrix_3x3_double *g, matrix_3x3_double *acc_prod);

inline extern void SU3_accumulate_right_product(matrix_3x3_double *acc_prod, const matrix_3x3_double *g);

inline extern void SU3_projection(matrix_3x3_double *x);

void SU3_decompose_algebra(const matrix_3x3_double *a, double *a_components);

#endif