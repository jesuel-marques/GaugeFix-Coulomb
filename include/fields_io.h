#ifndef CONFIGIO_H
#define CONFIGIO_H

#include <flags.h>
#include <types_io.h>
#include <lattice.h>

InGTMtrx  *getGaugetransfIn (InGTMtrx  * restrict G_in,  
                               const PosVec position);
OutGTMtrx *getGaugetransfOut(OutGTMtrx * restrict G_out, 
                               const PosVec position);

InCfgMtrx  *getLinkIn (InCfgMtrx  * restrict U_in,
                         const PosVec position, 
                         const LorentzIdx mu);
OutCfgMtrx *getLinkOut(OutCfgMtrx * restrict U_out,
                         const PosVec position,
                         const LorentzIdx mu);

int loadConfig (const Mtrx3x3 * restrict U,
                     char * config_filename),
    writeConfig(Mtrx3x3 * restrict U,
                     char * config_filename);

int loadGaugeTransf (Mtrx3x3 * restrict G,
                           char * gauge_transf_filename),
    writeGaugeTransf(Mtrx3x3 * restrict G,
                           char * gauge_transf_filename);

#endif  //CONFIGIO_H