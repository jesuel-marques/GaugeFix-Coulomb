#ifndef BICGSTAB_H
#define BICGSTAB_H

#include <types.h>

double BiCGStab(void (*operator)(Scalar *, Scalar *), Scalar *source, Scalar *inverse_column, double tolerance);

#endif  // BICGSTAB_H