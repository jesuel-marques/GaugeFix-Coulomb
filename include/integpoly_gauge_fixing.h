#ifndef INTEGPOLYGAUGEFIXING_H
#define INTEGPOLYGAUGEFIXING_H

#include <types.h>
#include <SU3_ops.h>

int integ_polyakov_gaugefix(Mtrx3x3 * restrict U, 
                            Mtrx3x3 * restrict G);

#endif  //INTEGPOLYGAUGEFIXING_H