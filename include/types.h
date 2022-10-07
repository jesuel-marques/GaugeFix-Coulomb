#ifndef TYPES_H
#define TYPES_H

#include <SU3_parameters.h>
#include <tgmath.h>

//  Geometric types definitions

typedef unsigned short pos_index;

typedef struct {
    pos_index i, j, k;
    pos_index t;
} pos_vec;  //	struct for position vectors


typedef unsigned short lorentz_idx;

typedef enum {REAR, FRONT} direction;

//  Data types definitions

typedef complex float  in_data_type;
typedef complex float  out_data_type;
typedef complex double work_data_type;
typedef work_data_type scalar;

// 2x2 types

typedef unsigned short SU2_color_idx;

typedef struct {
    scalar m[2 * 2];
} mtrx_2x2;

typedef struct {
    double m[4];
} mtrx_2x2_ck;

// 3x3 types

typedef unsigned short SU3_color_idx;
typedef unsigned short SU3_alg_idx;

typedef struct {
   complex double m[Nc * Nc];
} mtrx_3x3_double; 

typedef struct {
    complex float m[Nc * Nc];
} mtrx_3x3_float;

typedef struct {
   scalar m[Nc];
} color_3_vec;

typedef struct {
    scalar m[(Nc * Nc - 1) + 1];
} mtrx_SU3_alg;


typedef mtrx_3x3_double work_mtrx_data_type;
typedef work_mtrx_data_type mtrx_3x3;

typedef enum {R, S, T} submatrix;

#ifdef CONV_CFG_TO_WORKING_PRECISION
    typedef mtrx_3x3_float in_cfg_data_type;
#else
    typedef work_mtrx_data_type in_cfg_data_type;
#endif

#ifdef CONV_GT_TO_WORKING_PRECISION
    typedef mtrx_3x3_float  in_gt_data_type;
#else
    typedef mtrx_3x3_double in_gt_data_type;
#endif

#ifdef CONV_CFG_FROM_WORKING_PRECISION
    typedef mtrx_3x3_float  out_cfg_data_type;
#else
    typedef mtrx_3x3_double out_cfg_data_type;
#endif

#ifdef CONV_GT_FROM_WORKING_PRECISION
    typedef mtrx_3x3_float  out_gt_data_type;
#else
    typedef mtrx_3x3_double out_gt_data_type;
#endif

#endif