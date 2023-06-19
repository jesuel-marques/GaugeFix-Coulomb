#include <SU3_ops.h>
#include <fields.h>
#include <geometry.h>
#include <polyakov.h>

Mtrx3x3 polyakovLoop(PosVec position, Mtrx3x3* U) {
    Mtrx3x3 loop;
    t = 0;
    for (t = 0; t < Nt; t++) {
        position.pos[T_INDX] = t;
        getLink(U, position, T_INDX);
    }
}