INCLUDE_DIR := include
SOURCE_DIR := src
OBJECTS_DIR := obj

CC := gcc
CFLAGS := -I$(INCLUDE_DIR) -std=c99 -O4 -march=skylake-avx512 -mtune=skylake-avx512 -fopenmp -w

LIBS := -lm

DEPENDENCIES = $(wildcard $(INCLUDE_DIR)/*.h)

SOURCES := $(wildcard $(SOURCE_DIR)/*.c)

_OBJECTS = $(patsubst $(SOURCE_DIR)/%.c, %.o, $(SOURCES))
OBJECTS = $(patsubst %, $(OBJECTS_DIR)/%, $(_OBJECTS))


$(OBJECTS_DIR)/%.o: $(SOURCE_DIR)/%.c $(DEPENDENCIES) Makefile
	$(CC) -c -o $@ $< $(CFLAGS)

all: gauge_fix_coulomb plaquette_check

gauge_fix_coulomb: $(OBJECTS)
	$(CC) -o $@ $^ gauge_fix_coulomb.c $(CFLAGS) $(LIBS)

plaquette_check: $(OBJECTS)
	$(CC) -o $@ $^ plaquette_check.c $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(OBJECTS_DIR)/*.o *~ $(INCLUDE_DIR)/*~ $(SOURCE_DIR)/*~ $(OBJECTS_DIR)/*~ 