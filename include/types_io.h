#ifndef TYPES_IO_H
#define TYPES_IO_H

#include <flags.h>
#include <SU3_ops.h>
#include <tgmath.h>
#include <types.h>

/* Input scalar type */
typedef complex float  InScalar;
/* Output scalar type */
typedef complex float  OutScalar;

/* I'm aware of the rather complex structure of the input and output types below. */
/* The flags in flags.h decide which precision should be used for input and output */

#ifdef CONV_CFG_TO_WORKING_PRECISION
    /* Input configurations are in single precision */
    typedef Mtrx3x3Flt InCfgMtrx;
#else   //CONV_CFG_TO_WORKING_PRECISION
    /* Input configurations are in double precision */
    typedef Mtrx3x3Work InCfgMtrx;
#endif  //CONV_CFG_TO_WORKING_PRECISION

#ifdef CONV_GT_TO_WORKING_PRECISION
    /* Input gauge transformations are in single precision */    
    typedef Mtrx3x3Flt  InGTMtrx;
#else   //CONV_GT_TO_WORKING_PRECISION
    /* Input configurations are in double precision */
    typedef Mtrx3x3Dbl InGTMtrx;
#endif //CONV_GT_TO_WORKING_PRECISION

#ifdef CONV_CFG_FROM_WORKING_PRECISION
    /* Output configurations need to be in single precision */
    typedef Mtrx3x3Flt  OutCfgMtrx;
#else   //CONV_CFG_FROM_WORKING_PRECISION
    /* Output configurations need to be in double precision */
    typedef Mtrx3x3Dbl OutCfgMtrx;
#endif //CONV_CFG_FROM_WORKING_PRECISION

#ifdef CONV_GT_FROM_WORKING_PRECISION
    /* Output gauge-transformations need to be in single precision */
    typedef Mtrx3x3Flt  OutGTMtrx;
#else   //CONV_GT_FROM_WORKING_PRECISION
    /* Output gauge-transformations need to be in double precision */
    typedef Mtrx3x3Dbl OutGTMtrx;
#endif //CONV_GT_FROM_WORKING_PRECISION

#endif // TYPES_IO_H