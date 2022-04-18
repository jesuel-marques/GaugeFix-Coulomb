#ifndef SU3OPS_H
#define SU3OPS_H

typedef unsigned short SU3_color_index;
typedef unsigned short SU3_color_alg_index;

typedef struct su3_matrix{
    double complex m[3 * 3];
};

void SU3_print_matrix(const double complex *u, const char *name);

void SU3_copy(const double complex *u, double complex *u_copy);

void SU3_convert_fd(const float complex *u_float, double complex *u_double);

void SU3_convert_df(const double complex *u_double, float complex *u_float);

void SU3_set_to_null(double complex *u);

void SU3_set_to_identity(double complex *u);

void SU3_accumulate(const double complex *u, double complex *acc);

void SU3_subtraction(const double complex *u, const double complex *v, double complex *u_minus_v);

inline extern double complex SU3_trace(const double complex *u);

inline extern double complex SU3_determinant(const double complex *u);

inline extern void SU3_hermitean_conjugate(const double complex *u, double complex *u_dagger);

inline extern void SU3_multiplication_by_scalar(const double complex alpha, const double complex *u, double complex *alpha_times_u);
inline extern void SU3_substitution_multiplication_by_scalar(const double complex alpha, double complex *u);

inline extern void SU3_product(const double complex *u, const double complex *v, double complex *uv);

void SU3_product_three(const double complex *u, const double complex *v, const double complex *w, double complex *uvw);

void SU3_product_four(const double complex *u, const double complex *v, const double complex *w, const double complex *x, double complex *uvwx);

inline extern void SU3_accumulate_left_product(const double complex *g, double complex *acc_prod);

inline extern void SU3_accumulate_right_product(double complex *acc_prod, const double complex *g);

inline extern void SU3_projection(double complex *x);

void SU3_decompose_algebra(const double complex *a, double *a_components);

#endif