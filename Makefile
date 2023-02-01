INCLUDE_DIR := include
SOURCE_DIR := src
BUILD_DIR := build

RM = rm -f
CC := gcc
CFLAGS := -I$(INCLUDE_DIR) -std=c99 -O4 -march=skylake-avx512 -mtune=skylake-avx512 -fopenmp -w -DNUM_THREADS=`nproc`

WARNINGS:= -Wall -Wextra

LIBS := -lm


SOURCES := $(wildcard $(SOURCE_DIR)/*.c)

_OBJECTS = $(patsubst $(SOURCE_DIR)/%.c, %.o, $(SOURCES))
OBJECTS = $(patsubst %, $(BUILD_DIR)/%, $(_OBJECTS))

_DEPENDS = $(patsubst $(SOURCE_DIR)/%.c, %.d, $(SOURCES))
DEPENDS = $(patsubst %, $(BUILD_DIR)/%, $(_DEPENDS))

BINARIES = gauge_fix_coulomb

.PHONY: all clean

all: $(BINARIES)


gauge_fix_coulomb: $(OBJECTS) gauge_fix_coulomb.c
	$(CC) -o $@ $^ $(CFLAGS) $(WARNINGS) $(LIBS) 


gauge_fix_coulomb_mpi: $(OBJECTS) gauge_fix_coulomb_mpi.c
	$(CC) -o $@ $^ $(CFLAGS) $(WARNINGS) $(LIBS) 

-include $(DEPENDS)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.c Makefile
	$(CC) -c $< -o $@  $(CFLAGS) $(WARNINGS) -MMD -MP

clean:
	$(RM) $(BINARIES) $(OBJECTS) $(DEPENDS)