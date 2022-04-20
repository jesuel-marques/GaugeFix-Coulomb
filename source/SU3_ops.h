#ifndef SU3OPS_H
#define SU3OPS_H

#define elm(a,b)    (a) * Nc + (b)

typedef unsigned short SU3_color_index;
typedef unsigned short SU3_color_alg_index;

typedef struct {
    double complex m[Nc];
} color_3_vec;

typedef struct {
    double complex m[Nc * Nc];
} matrix_3x3_double;


typedef struct {
    float complex m[Nc * Nc];
} matrix_3x3_float;

typedef struct {
    double complex m[Nc * Nc - 1 + 1];
} matrix_SU3_alg;

void print_matrix_3x3(const matrix_3x3_double *u, const char *name);

void copy_3x3(const matrix_3x3_double *u, matrix_3x3_double *u_copy);

void convert_fd_3x3(const matrix_3x3_float *u_float, matrix_3x3_double *u_double),
     convert_df_3x3(const matrix_3x3_double *u_double, matrix_3x3_float *u_float);

inline void set_to_null_3x3(matrix_3x3_double *u),
            set_to_identity_3x3(matrix_3x3_double *u);

inline void accumulate_3x3(const matrix_3x3_double *u, matrix_3x3_double *acc);

void subtraction_3x3(const matrix_3x3_double *u, const matrix_3x3_double *v, matrix_3x3_double *u_minus_v);

inline extern double complex trace_3x3(const matrix_3x3_double *u);

inline extern double complex determinant_3x3(const matrix_3x3_double *u);

inline extern void SU3_hermitean_conjugate(const matrix_3x3_double *u, matrix_3x3_double *u_dagger);

inline extern void multiplication_by_scalar_3x3(const double complex alpha, const matrix_3x3_double *u, matrix_3x3_double *alpha_times_u);
inline extern void substitution_multiplication_by_scalar_3x3(const double complex alpha, matrix_3x3_double *u);

inline extern void product_3x3(const matrix_3x3_double *u, const matrix_3x3_double *v, matrix_3x3_double *uv);

void product_three_3x3(const matrix_3x3_double *u, const matrix_3x3_double *v, const matrix_3x3_double *w, matrix_3x3_double *uvw),
     product_four_3x3(const matrix_3x3_double *u, const matrix_3x3_double *v, const matrix_3x3_double *w, const matrix_3x3_double *x, matrix_3x3_double *uvwx);

inline extern void accumulate_left_product_3x3(const matrix_3x3_double *g, matrix_3x3_double *acc_prod),
                   accumulate_right_product_3x3(matrix_3x3_double *acc_prod, const matrix_3x3_double *g);

inline extern void projection_SU3(matrix_3x3_double *x);

void decompose_algebra_SU3(const matrix_3x3_double *a, matrix_SU3_alg *a_components);

#endif