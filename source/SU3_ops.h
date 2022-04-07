#ifndef SU3OPS_H
#define SU3OPS_H

void SU3_print_matrix(const float complex *u, const char *name);

void SU3_copy(const float complex *u, float complex *u_copy);

void SU3_set_to_null(float complex *u);

void SU3_set_to_identity(float complex *u);

void SU3_accumulate(const float complex *u, float complex *acc);

void SU3_subtraction(const float complex *u, const float complex *v, float complex *u_minus_v);

float complex SU3_trace(const float complex *u);

float complex SU3_determinant(const float complex *u);

void SU3_hermitean_conjugate(const float complex *u, float complex *u_dagger);

void SU3_multiplication_by_scalar(const float complex alpha, const float complex *u, float complex *alpha_times_u);
void SU3_substitution_multiplication_by_scalar(const float complex alpha, float complex *u);

void SU3_product(const float complex *u, const float complex *v, float complex *uv);

void SU3_product_three(const float complex *u, const float complex *v, const float complex *w, float complex *uvw);

void SU3_product_four(const float complex *u, const float complex *v, const float complex *w, const float complex *x, float complex *uvwx);

void SU3_accumulate_left_product(const float complex *g, float complex *acc_prod);

void SU3_accumulate_right_product(float complex *acc_prod, const float complex *g);

void SU3_projection(float complex *x);

void SU3_decompose_algebra(const float complex *a, float *a_components);

#endif