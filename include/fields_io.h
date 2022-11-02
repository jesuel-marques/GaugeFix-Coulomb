#ifndef CONFIGIO_H
#define CONFIGIO_H

#include <flags.h>
#include <io_types.h>
#include <types.h>

InGTMtrx  *get_gaugetransf_in (InGTMtrx  * restrict G_in,  
                               const PosVec position);
OutGTMtrx *get_gaugetransf_out(OutGTMtrx * restrict G_out, 
                               const PosVec position);

InCfgMtrx  *get_link_in (InCfgMtrx  * restrict U_in,
                         const PosVec position, 
                         const LorentzIdx mu);
OutCfgMtrx *get_link_out(OutCfgMtrx * restrict U_out,
                         const PosVec position,
                         const LorentzIdx mu);

int SU3_load_config (Mtrx3x3 * restrict U,
                     char * config_filename),
    SU3_write_config(Mtrx3x3 * restrict U,
                     char * config_filename);

int SU3_load_gauge_transf (Mtrx3x3 * restrict G,
                           char * gauge_transf_filename),
    SU3_write_gauge_transf(Mtrx3x3 * restrict G,
                           char * gauge_transf_filename);

#endif  //CONFIGIO_H