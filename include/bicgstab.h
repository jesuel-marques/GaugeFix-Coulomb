#ifndef BICGSTAB_H
#define BICGSTAB_H

#include <stddef.h>
#include <types.h>

double BiCGStab(void (*operator)(Scalar *, Scalar *), Scalar *source,
                Scalar *inverse_column, double tolerance, size_t sizeof_vector);

#endif  // BICGSTAB_H