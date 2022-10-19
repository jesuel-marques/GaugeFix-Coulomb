#ifndef GAUGEFIXING_H
#define GAUGEFIXING_H

#include <types.h>

void SU3_global_update_U(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G);

void SU3_update_sub_LosAlamos(mtrx_3x3 * restrict w, submatrix sub);

void SU3_LosAlamos_common_block(mtrx_3x3 * restrict w, 
                                mtrx_3x3 * restrict total_update);

double SU3_calculate_F(mtrx_3x3 * restrict U);

double SU3_calculate_theta(mtrx_3x3 * restrict U);

double SU3_calculate_e2(mtrx_3x3 * restrict U);

int SU3_gauge_fix(mtrx_3x3 * restrict U,  mtrx_3x3 * restrict G, const unsigned short config);

void init_gauge_transformation(mtrx_3x3 * restrict G);

#endif