#include <convert_precision_io.h>
#include <fields.h>
#include <fields_io.h>
#include <flags.h>
#include <io_types.h>
#include <lattice.h>
#include <misc.h>
#include <SU3_ops.h>
#include <types.h>


extern geometric_parameters lattice_param;

#ifdef CONV_CFG_TO_WORKING_PRECISION
static void convert_in_work_cfg_3x3(const InCfgMtrx * restrict u_in, 
                                          WorkMtrx  * restrict u_work) {
    MtrxIdx3 a, b;
    LOOP_3X3(a, b){

        u_work -> m[ELEM_3X3(a, b)] = (WorkScalarType) u_in -> m[ELEM_3X3(a, b)];

    }
}
#endif  //CONV_CFG_TO_WORKING_PRECISION


#ifdef CONV_CFG_FROM_WORKING_PRECISION
static void convert_work_out_cfg_3x3(const WorkMtrx   * restrict u_work, 
                                           OutCfgMtrx * restrict u_out) {
    MtrxIdx3 a, b;
    LOOP_3X3(a, b){

        u_out -> m[ELEM_3X3(a, b)] = (OutScalar) u_work -> m[ELEM_3X3(a, b)];

    }
}
#endif  //CONV_CFG_FROM_WORKING_PRECISION


#ifdef CONV_GT_TO_WORKING_PRECISION
static void convert_in_work_gt_3x3(const InGTMtrx * restrict g_in, 
                                         WorkMtrx * restrict g_work) {
    MtrxIdx3 a, b;
    LOOP_3X3(a, b){

        g_work -> m[ELEM_3X3(a, b)] = (WorkScalarType) g_in -> m[ELEM_3X3(a, b)];

    }
}
#endif  //CONV_GT_TO_WORKING_PRECISION


#ifdef CONV_GT_FROM_WORKING_PRECISION
static void convert_work_out_gt_3x3(const WorkMtrx  * restrict g_work, 
                                          OutGTMtrx * restrict g_out) {
    MtrxIdx3 a, b;
    LOOP_3X3(a, b){

        g_out -> m[ELEM_3X3(a, b)] = (OutScalar) g_work -> m[ELEM_3X3(a, b)];

    }
}
#endif  // CONV_GT_FROM_WORKING_PRECISION


#ifdef CONV_CFG_TO_WORKING_PRECISION
void SU3_convert_cfg_in_work(InCfgMtrx * restrict U_in, 
                              WorkMtrx * restrict U_work) {
    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t){
        PosVec position;
        position.t = t;
        LOOP_SPATIAL(position){
            LorentzIdx mu;
            LOOP_LORENTZ(mu){

                convert_in_work_cfg_3x3(get_link_in(U_in  , position, mu), 
                                        get_link   (U_work, position, mu));

            }
        }
    }
}
#endif  //CONV_CFG_TO_WORKING_PRECISION


#ifdef CONV_CFG_FROM_WORKING_PRECISION
void SU3_convert_cfg_work_out(WorkMtrx   * restrict U_work, 
                              OutCfgMtrx * restrict U_out) {
    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t){
        PosVec position;
        position.t = t;
        LOOP_SPATIAL(position){
            LorentzIdx mu;
            LOOP_LORENTZ(mu) {

                convert_work_out_cfg_3x3(get_link    (U_work, position, mu), 
                                         get_link_out(U_out,  position, mu));

            }
        }
    }
}
#endif  //CONV_CFG_FROM_WORKING_PRECISION


#ifdef CONV_GT_TO_WORKING_PRECISION
 void SU3_convert_gt_in_work(InGTMtrx * restrict G_in, 
                             WorkMtrx * restrict G_work) {
    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t){
        PosVec position;
        position.t = t;
        LOOP_SPATIAL(position){

            convert_in_work_gt_3x3(get_gaugetransf_in(G_in,   position), 
                                   get_gaugetransf   (G_work, position));                    

        }
    }
}
#endif  //CONV_GT_TO_WORKING_PRECISION


#ifdef CONV_GT_FROM_WORKING_PRECISION
void SU3_convert_gt_work_out(WorkMtrx  *G_work, 
                             OutGTMtrx *G_out) {
    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t){
        PosVec position;
        position.t = t;
        LOOP_SPATIAL(position){

            convert_work_out_gt_3x3(get_gaugetransf    (G_work, position), 
                                    get_gaugetransf_out(G_out,  position));                       

        }
    }
}
#endif  //CONV_GT_FROM_WORKING_PRECISION