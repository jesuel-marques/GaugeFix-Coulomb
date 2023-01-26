#ifndef CONVERT_PRECISION_IO_H
#define  CONVERT_PRECISION_IO_H


#include <flags.h>
#include <types_io.h>
#include <lattice.h>
#include <types.h>
#include <SU3_ops.h>


#ifdef CONV_CFG_TO_WORKING_PRECISION

void convertCfg_in_work(InCfgMtrx * restrict U_in, Mtrx3x3Work  * restrict U_work);

#endif  //CONV_CFG_TO_WORKING_PRECISION

#ifdef CONV_GT_TO_WORKING_PRECISION

void convertGT_in_work(InGTMtrx * restrict G_in, Mtrx3x3Work * restrict G_work);

#endif  //CONV_GT_TO_WORKING_PRECISION


#ifdef CONV_CFG_FROM_WORKING_PRECISION

void convertCfg_work_out(Mtrx3x3Work * restrict U_work, OutCfgMtrx * restrict U_out);

#endif  //CONV_CFG_FROM_WORKING_PRECISION

#ifdef CONV_GT_FROM_WORKING_PRECISION

void convertGT_work_out(Mtrx3x3Work *G_work, OutGTMtrx *G_out);

#endif  //CONV_GT_FROM_WORKING_PRECISION

#endif  //CONVERT_PRECISION_IO_H