#include <convert_precision_io.h>
#include <fields.h>
#include <fields_io.h>
#include <flags.h>
#include <types_io.h>
#include <lattice.h>
#include <misc.h>
#include <SU3_ops.h>
#include <types.h>

#ifdef CONV_CFG_TO_WORKING_PRECISION
/* Converts a SU(3) gauge link from the input precision to the internal working 
   precision. */
static void convertCfg_in_work_3x3(const InCfgMtrx * restrict u_in, 
                                         WorkMtrx  * restrict u_work) {

    /*
	 * Calls:
	 * =====
	 *
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * InCfgMtrx * u_in:    SU(3) gauge-link in input precision,
     * WorkMtrx  * u_work:  SU(3) gauge-link in internal working precision.
     * 
	 * Returns:
	 * =======
	 * 
     */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {

        u_work -> m[ELEM_3X3(a, b)] = (WorkScalarType) u_in -> m[ELEM_3X3(a, b)];

    }
}
#endif  //CONV_CFG_TO_WORKING_PRECISION


#ifdef CONV_CFG_FROM_WORKING_PRECISION
/* Converts a SU(3) gauge link from the internal working precision to the output 
   precision. */
static void convertCfg_work_out_3x3(const WorkMtrx * restrict u_work, 
                                    OutCfgMtrx * restrict u_out) {

    /*
	 * Calls:
	 * =====
	 *
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * WorkMtrx * u_work:   SU(3) gauge-link in internal working precision,
     * OutCfgMtrx * u_out:  SU(3) gauge-link in output precision.
     * 
	 * Returns:
	 * =======
	 * 
     */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        u_out -> m[ELEM_3X3(a, b)] = (OutScalar) u_work -> m[ELEM_3X3(a, b)];
    }
}
#endif  //CONV_CFG_FROM_WORKING_PRECISION


#ifdef CONV_GT_TO_WORKING_PRECISION
/* Converts a SU(3) gauge-transformation matrix from the input precision to the
   internal working precision. */
static void convertGT_in_work_3x3(const InGTMtrx * restrict g_in, 
                                        WorkMtrx * restrict g_work) {

    /* 
	 * Calls:
	 * =====
	 *
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * InGTMtrx * g_in:     SU(3) gauge-transformation matrix in input precision,
     * WorkMtrx  * g_work:  SU(3) gauge-transformation matrix in working precision.
     * 
	 * Returns:
	 * =======
	 * 
     */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        g_work -> m[ELEM_3X3(a, b)] = (WorkScalarType) g_in -> m[ELEM_3X3(a, b)];
    }
}
#endif  //CONV_GT_TO_WORKING_PRECISION


#ifdef CONV_GT_FROM_WORKING_PRECISION
/* Converts a SU(3) gauge-transformation matrix from the internal working precision 
   to the output precision. */
static void convertGT_work_out_3x3(const WorkMtrx  * restrict g_work, 
                                         OutGTMtrx * restrict g_out) {

    /* 
	 * Calls:
	 * =====
	 *
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * WorkMtrx * g_work:   SU(3) gauge-transformation matrix in working precision,
     * OutGTMtrx * g_out:   SU(3) gauge-transformation matrix in output precision.
     * 
	 * Returns:
	 * =======
	 * 
     */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        g_out -> m[ELEM_3X3(a, b)] = (OutScalar) g_work -> m[ELEM_3X3(a, b)];
    }
}
#endif  // CONV_GT_FROM_WORKING_PRECISION


#ifdef CONV_CFG_TO_WORKING_PRECISION
/* Converts a SU(3) gauge field from the input precision to the internal working 
   precision. */
void convertCfg_in_work(InCfgMtrx * restrict U_in, 
                         WorkMtrx * restrict U_work) {

    /* 
	 * Calls:
	 * =====
     * getLink, getLinkIn,
	 * convertCfg_in_work_3x3.
     * 
     * Macros:
	 * ======
     * LOOP_TEMPORAL_PARALLEL, T_INDX, LOOP_SPATIAL, LOOP_LORENTZ.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * InCfgMtrx * U_in:    SU(3) gauge field in the input precision,
     * WorkMtrx * U_work:   SU(3) gauge field in the interal working precision.
     * 
	 * Returns:
	 * =======
	 * 
     */

    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t) {
        PosVec position;
        position.pos[T_INDX] = t;
        LOOP_SPATIAL(position) {
            LorentzIdx mu;
            LOOP_LORENTZ(mu) {
                convertCfg_in_work_3x3(getLinkIn(U_in  , position, mu), 
                                       getLink  (U_work, position, mu));
            }
        }
    }
}
#endif  //CONV_CFG_TO_WORKING_PRECISION


#ifdef CONV_CFG_FROM_WORKING_PRECISION
/* Converts a SU(3) gluon field from the internal working precision to  output 
   precision. */
void convertCfg_work_out(WorkMtrx * restrict U_work, OutCfgMtrx * restrict U_out) {
    
    /* 
	 * Calls:
	 * =====
     * getLink, getLinkOut,
	 * convertCfg_work_out_3x3.
     * 
     * Macros:
	 * ======
     * LOOP_TEMPORAL_PARALLEL, LOOP_SPATIAL, LOOP_LORENTZ.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * WorkMtrx * U_work:   SU(3) gluon field in the internal working precision,
     * OutCfgMtrx * U_out:  SU(3) gluon field in the output precision.
     * 
	 * Returns:
	 * =======
	 * 
     */

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
/* Converts a SU(3) gauge-transformation field from the input precision to the
   internal working precision.  */
 void convertGT_in_work(InGTMtrx * restrict G_in, WorkMtrx * restrict G_work) {
    
    /* 
	 * Calls:
	 * =====
     * getGaugetransf, getGaugetransfIn,
	 * convertGT_in_work_3x3.
     *
     * Macros:
	 * ======
     * LOOP_TEMPORAL_PARALLEL, LOOP_SPATIAL.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * InGTMtrx * G_in:     SU(3) gauge-transformation field in input precision,
     * WorkMtrx * G_work:   SU(3) gauge-transformation field in internal working 
     *                      precision.
     * 
	 * Returns:
	 * =======
	 * 
     */

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
/* Converts a SU(3) gauge-transformation field from the internal working precision 
   to the output precision. */
void convertGT_work_out(WorkMtrx  *G_work, 
                        OutGTMtrx *G_out) {

    /*  
	 * Calls:
	 * =====
     * getGaugetransf, getGaugetransfOut,
	 * convertGT_work_out_3x3.
     *
     * Macros:
	 * ======
     * LOOP_TEMPORAL_PARALLEL, LOOP_SPATIAL.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * WorkMtrx * G_work:   SU(3) gauge-transformation field in internal working 
     *                      precision,
     * OutGTMtrx * G_out:   SU(3) gauge-transformation in output precision.
     * 
	 * Returns:
	 * =======
	 * 
     */

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