#ifndef IO_TYPES_H
#define IO_TYPES_H

#include <flags.h>
#include <SU3_ops.h>
#include <tgmath.h>
#include <types.h>

typedef complex float  InScalar;
typedef complex float  OutScalar;

#ifdef CONV_CFG_TO_WORKING_PRECISION
    typedef Mtrx3x3Flt InCfgMtrx;
#else   //CONV_CFG_TO_WORKING_PRECISION
    typedef WorkMtrx InCfgMtrx;
#endif  //CONV_CFG_TO_WORKING_PRECISION

#ifdef CONV_GT_TO_WORKING_PRECISION
    typedef Mtrx3x3Flt  InGTMtrx;
#else   //CONV_GT_TO_WORKING_PRECISION
    typedef Mtrx3x3Dbl InGTMtrx;
#endif //CONV_GT_TO_WORKING_PRECISION

#ifdef CONV_CFG_FROM_WORKING_PRECISION
    typedef Mtrx3x3Flt  OutCfgMtrx;
#else   //CONV_CFG_FROM_WORKING_PRECISION
    typedef Mtrx3x3Dbl OutCfgMtrx;
#endif //CONV_CFG_FROM_WORKING_PRECISION

#ifdef CONV_GT_FROM_WORKING_PRECISION
    typedef Mtrx3x3Flt  OutGTMtrx;
#else   //CONV_GT_FROM_WORKING_PRECISION
    typedef Mtrx3x3Dbl OutGTMtrx;
#endif //CONV_GT_FROM_WORKING_PRECISION

#endif // IO_TYPES_H