#ifndef SU3OPS_H
#define SU3OPS_H

#define elm(a,b)    (a) * Nc + (b)

typedef unsigned short SU3_color_idx;
typedef unsigned short SU3_alg_idx;

typedef struct {
    double complex m[Nc];
} color_3_vec;

typedef struct {
    double complex m[Nc * Nc];
} mtrx_3x3_double; 


typedef struct {
    float complex m[Nc * Nc];
} mtrx_3x3_float;

typedef struct {
    double complex m[Nc * Nc - 1 + 1];
} matrix_SU3_alg;

void print_matrix_3x3(const mtrx_3x3_double *u, const char *name);

void copy_3x3(const mtrx_3x3_double *u, mtrx_3x3_double *u_copy);

void convert_fd_3x3(const mtrx_3x3_float *u_float, mtrx_3x3_double *u_double),
     convert_df_3x3(const mtrx_3x3_double *u_double, mtrx_3x3_float *u_float);

void set_null_3x3(mtrx_3x3_double *u),
     set_identity_3x3(mtrx_3x3_double *u);

void accumulate_3x3(const mtrx_3x3_double *u, mtrx_3x3_double *acc);

void subtraction_3x3(const mtrx_3x3_double *u, const mtrx_3x3_double *v, mtrx_3x3_double *u_minus_v);

double complex trace_3x3(const mtrx_3x3_double *u),
         determinant_3x3(const mtrx_3x3_double *u);

void SU3_hermitean_conjugate(const mtrx_3x3_double *u, mtrx_3x3_double *u_dagger);

void multiplication_by_scalar_3x3(const double complex alpha, const mtrx_3x3_double *u, mtrx_3x3_double *alpha_times_u);
void substitution_multiplication_by_scalar_3x3(const double complex alpha, mtrx_3x3_double *u);

void product_3x3(const mtrx_3x3_double *u, const mtrx_3x3_double *v, mtrx_3x3_double *uv);

void product_three_3x3(const mtrx_3x3_double *u, const mtrx_3x3_double *v, const mtrx_3x3_double *w, mtrx_3x3_double *uvw),
     product_four_3x3(const mtrx_3x3_double *u, const mtrx_3x3_double *v, const mtrx_3x3_double *w, const mtrx_3x3_double *x, mtrx_3x3_double *uvwx);

void accumulate_left_prod_3x3(const mtrx_3x3_double *g, mtrx_3x3_double *acc_prod),
     accumulate_right_prod_3x3(mtrx_3x3_double *acc_prod, const mtrx_3x3_double *g);

void projection_SU3(mtrx_3x3_double *x);

void decompose_algebra_SU3(const mtrx_3x3_double *a, matrix_SU3_alg *a_components);

#endif