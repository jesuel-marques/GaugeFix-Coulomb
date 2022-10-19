INCLUDE_DIR := include
SOURCE_DIR := src
BUILD_DIR := build

LAPACK_ROOT_DIR := /home/postgrad/jesuel/lapack/lapack-3.10.1
LAPACK_INCLUDE_DIR := $(LAPACK_ROOT_DIR)/LAPACKE/include
LAPACK_OBJS := $(LAPACK_ROOT_DIR)/liblapacke.a $(LAPACK_ROOT_DIR)/liblapack.a $(LAPACK_ROOT_DIR)/librefblas.a $(LAPACK_ROOT_DIR)/libtmglib.a

CC := mpiicc
CFLAGS := -I$(INCLUDE_DIR) -I$(LAPACK_INCLUDE_DIR) -std=c99 -O3 -ipo -xHASWELL -axSKYLAKE,CASCADELAKE,TIGERLAKE -qopt-zmm-usage=high -qopenmp -DMPI_CODE

WARNINGS:= -Wall -Wextra

LIBS := -lm -lgfortran

SOURCES := $(wildcard $(SOURCE_DIR)/*.c)

_OBJECTS = $(patsubst $(SOURCE_DIR)/%.c, %.o, $(SOURCES))
OBJECTS = $(patsubst %, $(BUILD_DIR)/%, $(_OBJECTS))

_DEPENDS = $(patsubst $(SOURCE_DIR)/%.c, %.d, $(SOURCES))
DEPENDS = $(patsubst %, $(BUILD_DIR)/%, $(_DEPENDS))

BINARIES = gauge_fix_coulomb plaquette_check test_matrix_power

.PHONY: all clean

all: $(BINARIES)

gauge_fix_coulomb: $(OBJECTS)
	$(CC) -o $@ $@.c $^ $(LAPACK_OBJS) $(CFLAGS) $(WARNINGS) $(LIBS) 

plaquette_check: $(OBJECTS)
	$(CC) -o $@ $@.c $^ $(LAPACK_OBJS) $(CFLAGS) $(WARNINGS) $(LIBS)

test_matrix_power: $(OBJECTS)
	$(CC) -o $@ $@.c $^ $(LAPACK_OBJS) $(CFLAGS) $(WARNINGS) $(LIBS)

-include $(DEPENDS)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.c Makefile
	$(CC) -c $< -o $@  $(CFLAGS) $(WARNINGS) -MMD -MP

clean:
	rm -f $(BINARIES) $(BUILD_DIR)/*.[od] *~ */*~ 