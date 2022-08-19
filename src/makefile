IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -std=c99 -O4 -march=skylake -fopenmp -w

ODIR=obj

LIBS=-lm

_DEPS = lattice.h SU2_ops.h SU3_ops.h gauge_fixing.h fourvector_field.h config_io.h integpoly_gauge_fixing.h math_ops.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = gauge_fix_coulomb.o lattice.o SU2_ops.o SU3_ops.o gauge_fixing.o fourvector_field.o config_io.o integpoly_gauge_fixing.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gauge_fix_coulomb: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 