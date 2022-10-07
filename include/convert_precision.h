#ifndef CONVERT_PRECISION_H
#define  CONVERT_PRECISION_H

#include <types.h>

#ifdef CONV_CFG_TO_WORKING_PRECISION
void SU3_convert_cfg_in_work(in_cfg_data_type * restrict U_in, work_mtrx_data_type * restrict U_work) ;

#endif

#ifdef CONV_GT_TO_WORKING_PRECISION

void SU3_convert_gt_in_work(in_gt_data_type * restrict G_in, work_mtrx_data_type * restrict G_work);
#endif


#ifdef CONV_CFG_FROM_WORKING_PRECISION

void SU3_convert_cfg_work_out(work_mtrx_data_type * restrict U_work, out_cfg_data_type * restrict U_out) ;
#endif

#ifdef CONV_GT_FROM_WORKING_PRECISION

void SU3_convert_gt_work_out(work_mtrx_data_type *G_work, out_gt_data_type *G_out);

#endif

#endif