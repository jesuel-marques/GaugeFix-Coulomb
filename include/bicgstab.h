#ifndef BICGSTAB_H
#define BICGSTAB_H

#include <stddef.h>
#include <types.h>

double BiCGStab(Scalar *(*operator)(Scalar *, Scalar *), Scalar *source,
                Scalar *inverse_column, double tolerance, size_t sizeof_vector);

double CG(Scalar *(*operator)(Scalar *, Scalar *), Scalar *restrict source, Scalar *restrict inverse_column, const double tolerance, const size_t sizeof_vector);

#endif  // BICGSTAB_H