#include <convert_precision_io.h>
#include <fields.h>
#include <fields_io.h>
#include <flags.h>
#include <io_types.h>
#include <lattice.h>
#include <misc.h>
#include <SU3_ops.h>
#include <types.h>

#ifdef CONV_CFG_TO_WORKING_PRECISION
static void convertCfg_in_work_3x3(const InCfgMtrx * restrict u_in, 
                                         WorkMtrx  * restrict u_work) {
    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {

        u_work -> m[ELEM_3X3(a, b)] = (WorkScalarType) u_in -> m[ELEM_3X3(a, b)];

    }
}
#endif  //CONV_CFG_TO_WORKING_PRECISION


#ifdef CONV_CFG_FROM_WORKING_PRECISION
static void convertCfg_work_out_3x3(const WorkMtrx   * restrict u_work, 
                                          OutCfgMtrx * restrict u_out) {
    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {

        u_out -> m[ELEM_3X3(a, b)] = (OutScalar) u_work -> m[ELEM_3X3(a, b)];

    }
}
#endif  //CONV_CFG_FROM_WORKING_PRECISION


#ifdef CONV_GT_TO_WORKING_PRECISION
static void convertGT_in_work_3x3(const InGTMtrx * restrict g_in, 
                                        WorkMtrx * restrict g_work) {
    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {

        g_work -> m[ELEM_3X3(a, b)] = (WorkScalarType) g_in -> m[ELEM_3X3(a, b)];

    }
}
#endif  //CONV_GT_TO_WORKING_PRECISION


#ifdef CONV_GT_FROM_WORKING_PRECISION
static void convertGT_work_out_3x3(const WorkMtrx  * restrict g_work, 
                                         OutGTMtrx * restrict g_out) {
    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {

        g_out -> m[ELEM_3X3(a, b)] = (OutScalar) g_work -> m[ELEM_3X3(a, b)];

    }
}
#endif  // CONV_GT_FROM_WORKING_PRECISION


#ifdef CONV_CFG_TO_WORKING_PRECISION
void convertCfg_in_work(InCfgMtrx * restrict U_in, 
                         WorkMtrx * restrict U_work) {
    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t) {
        PosVec position;
        position.t = t;
        LOOP_SPATIAL(position) {
            LorentzIdx mu;
            LOOP_LORENTZ(mu) {

                convertCfg_in_work_3x3(getLinkIn(U_in  , position, mu), 
                                        getLink (U_work, position, mu));

            }
        }
    }
}
#endif  //CONV_CFG_TO_WORKING_PRECISION


#ifdef CONV_CFG_FROM_WORKING_PRECISION
void convertCfg_work_out(WorkMtrx   * restrict U_work, 
                          OutCfgMtrx * restrict U_out  ) {
    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t) {
        PosVec position;
        position.t = t;
        LOOP_SPATIAL(position) {
            LorentzIdx mu;
            LOOP_LORENTZ(mu) {

                convertCfg_work_out_3x3(getLink    (U_work, position, mu), 
                                        getLinkOut (U_out,  position, mu));

            }
        }
    }
}
#endif  //CONV_CFG_FROM_WORKING_PRECISION


#ifdef CONV_GT_TO_WORKING_PRECISION
 void convertGT_in_work(InGTMtrx * restrict G_in, 
                         WorkMtrx * restrict G_work) {
    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t) {
        PosVec position;
        position.t = t;
        LOOP_SPATIAL(position) {

            convertGT_in_work_3x3(getGaugetransfIn(G_in,   position), 
                                  getGaugetransf  (G_work, position));                    

        }
    }
}
#endif  //CONV_GT_TO_WORKING_PRECISION


#ifdef CONV_GT_FROM_WORKING_PRECISION
void convertGT_work_out(WorkMtrx  *G_work, 
                         OutGTMtrx *G_out) {
    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t) {
        PosVec position;
        position.t = t;
        LOOP_SPATIAL(position) {

            convertGT_work_out_3x3(getGaugetransf    (G_work, position), 
                                   getGaugetransfOut (G_out,  position));                       

        }
    }
}
#endif  //CONV_GT_FROM_WORKING_PRECISION