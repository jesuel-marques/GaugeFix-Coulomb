#ifndef CONVERT_PRECISION_H
#define  CONVERT_PRECISION_H


#include <flags.h>
#include <io_types.h>
#include <types.h>


#ifdef CONV_CFG_TO_WORKING_PRECISION
void SU3_convert_cfg_in_work(InCfgMtrx * restrict U_in, 
                             WorkMtrx  * restrict U_work);

#endif  //CONV_CFG_TO_WORKING_PRECISION

#ifdef CONV_GT_TO_WORKING_PRECISION

void SU3_convert_gt_in_work(InGTMtrx * restrict G_in, 
                            WorkMtrx * restrict G_work);
#endif  //CONV_GT_TO_WORKING_PRECISION


#ifdef CONV_CFG_FROM_WORKING_PRECISION

void SU3_convert_cfg_work_out(WorkMtrx   * restrict U_work,
                              OutCfgMtrx * restrict U_out );
#endif  //CONV_CFG_FROM_WORKING_PRECISION

#ifdef CONV_GT_FROM_WORKING_PRECISION

void SU3_convert_gt_work_out(WorkMtrx  *G_work,
                             OutGTMtrx *G_out );

#endif  //CONV_GT_FROM_WORKING_PRECISION

#endif  //CONVERT_PRECISION_H