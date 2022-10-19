#ifndef INTEGPOLYGAUGEFIXING_H
#define INTEGPOLYGAUGEFIXING_H

#include <types.h>

int integ_polyakov_gauge_fix(mtrx_3x3 * restrict U, 
                            mtrx_3x3 * restrict G, const unsigned short config_nr);

#endif