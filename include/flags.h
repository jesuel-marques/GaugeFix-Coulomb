#ifndef FLAGS_H
#define FLAGS_H

/* Uncomment if conversion from big to little endian is needed in config being read */
#define     NEED_BYTE_SWAP_IN                   

/* Uncomment if conversion from big to little endian is needed in config being written*/
// #define  NEED_BYTE_SWAP_OUT                  


/* Uncomment if reading config in float and need to convert to double precision */
#define     CONV_CFG_TO_WORKING_PRECISION   

/*  Uncomment if writing config in float and need to convert from double precision */
// #define     CONV_CFG_FROM_WORKING_PRECISION


/*  Uncomment if reading gauge-transf in float and 
    need to convert to double precision */
// #define  CONV_GT_TO_WORKING_PRECISION

/*  Uncomment if writing gauge-transf in float and 
    need to convert from double precision */
// #define  CONV_GT_FROM_WORKING_PRECISION


// #define CHECK_POSITION_BOUNDS
//  Uncomment to check that only positions within bounds are being accessed

#endif  //FLAGS_H