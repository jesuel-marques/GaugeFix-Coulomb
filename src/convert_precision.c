#include <types.h>
#include <fields.h>
#include <convert_precision.h>
#include <SU3_ops.h>
#include <SU3_parameters.h>
#include <misc.h>

static void convert_in_work_cfg_3x3(const in_cfg_data_type * restrict u_in, 
                                work_mtrx_data_type * restrict u_work) {

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            u_work -> m[ELM3x3(a, b)] = (work_data_type) u_in -> m[ELM3x3(a, b)];

        }
    }
}

static void convert_work_out_cfg_3x3(const work_mtrx_data_type * restrict u_work, 
                                      out_cfg_data_type * restrict u_out) {

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            u_out -> m[ELM3x3(a, b)] = (out_data_type) u_work -> m[ELM3x3(a, b)];

        }
    }
}


static void convert_in_work_gt_3x3(const in_gt_data_type * restrict g_in, 
                              work_mtrx_data_type * restrict g_work) {

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            g_work -> m[ELM3x3(a, b)] = (work_data_type) g_in -> m[ELM3x3(a, b)];

        }
    }
}

static void convert_work_out_gt_3x3(const work_mtrx_data_type * restrict g_work, 
                                      out_gt_data_type * restrict g_out) {

    for (SU3_color_idx  a = 0; a < Nc; a++) {
        for (SU3_color_idx  b = 0; b < Nc; b++) {

            g_out -> m[ELM3x3(a, b)] = (out_data_type) g_work -> m[ELM3x3(a, b)];

        }
    }
}


#ifdef CONV_CFG_TO_WORKING_PRECISION
void SU3_convert_cfg_in_work(in_cfg_data_type * restrict U_in, work_mtrx_data_type * restrict U_work) {
    OMP_PARALLEL_FOR
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {
                            convert_in_work_cfg_3x3(get_link_in(U_in, position, mu), get_link(U_work, position, mu));
                        }
                    }
                }
            }
        }
}

#endif

#ifdef CONV_GT_TO_WORKING_PRECISION

 void SU3_convert_gt_in_work(in_gt_data_type * restrict G_in, work_mtrx_data_type * restrict G_work) {
        OMP_PARALLEL_FOR
            // Paralelizing by slicing the time extent
            for (pos_index t = 0; t < N_T; t++) {
                pos_vec position;
                position.t = t;
                for (position.k = 0; position.k < N_SPC; position.k++) {
                    for (position.j = 0; position.j < N_SPC; position.j++) {
                        for (position.i = 0; position.i < N_SPC; position.i++) {                    
                                convert_in_work_gt_3x3(get_gaugetransf_in(G_in, position), get_gaugetransf(G_work, position));                    
                        }
                    }
                }
            }
    }
#endif


#ifdef CONV_CFG_FROM_WORKING_PRECISION

void SU3_convert_cfg_work_out(work_mtrx_data_type * restrict U_work, out_cfg_data_type * restrict U_out) {
    OMP_PARALLEL_FOR
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {
                            convert_work_out_cfg_3x3(get_link(U_work, position, mu), get_link_out(U_out, position, mu));
                        }
                    }
                }
            }
        }
}

#endif

#ifdef CONV_GT_FROM_WORKING_PRECISION

void SU3_convert_gt_work_out(work_mtrx_data_type *G_work, out_gt_data_type *G_out) {
    OMP_PARALLEL_FOR
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {                    
                            convert_work_out_gt_3x3(get_gaugetransf(G_work, position), get_gaugetransf_out(G_out, position));                    
                    }
                }
            }
        }
}

#endif
