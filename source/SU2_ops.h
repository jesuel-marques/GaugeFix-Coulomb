#ifndef SU2OPS_H
#define SU2OPS_H

void SU2_copy(const float* u, float* u_copy);

void SU2_set_to_null(float complex* u);
void SU2_set_to_identity(float complex* u);

void SU2_accumulate(const float* u, float* acc);

void SU2_subtraction(const float* u, const float* v, float* u_minus_v);

float SU2_determinant(const float* u);

void SU2_hermitean_conjugate(const float* u, float* u_dagger);

void SU2_multiplication_by_scalar(const float* u, const float alpha, float* alpha_times_u);

void SU2_product(const float* u, const float* v, float* uv);

void SU2_product_three(const float* u, const float* v, const float* w, float* uvw);

void SU2_product_four(const float* u, const float* v, const float* w, const float* x, float* uvwx);

void SU2_projection(float* u);

#endif