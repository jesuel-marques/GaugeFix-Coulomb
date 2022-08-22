INCLUDE_DIR := include
SOURCE_DIR := src
BUILD_DIR := build
DEPENDS_DIR := dep

CC := gcc
CFLAGS := -I$(INCLUDE_DIR) -std=c99 -O4 -march=skylake-avx512 -mtune=skylake-avx512 -fopenmp -w

WARNINGS:= -Wall -Wextra

LIBS := -lm

DEPENDENCIES = $(wildcard $(INCLUDE_DIR)/*.h)

SOURCES := $(wildcard $(SOURCE_DIR)/*.c)

_OBJECTS = $(patsubst $(SOURCE_DIR)/%.c, %.o, $(SOURCES))
OBJECTS = $(patsubst %, $(BUILD_DIR)/%, $(_OBJECTS))

_DEPENDS = $(patsubst $(SOURCE_DIR)/%.c, %.d, $(SOURCES))
DEPENDS = $(patsubst %, $(BUILD_DIR)/%, $(_DEPENDS))

BINARIES = gauge_fix_coulomb plaquette_check

.PHONY: all clean

all: $(BINARIES)

gauge_fix_coulomb: $(OBJECTS)
	$(CC) -o $@ $@.c $^ $(CFLAGS) $(WARNINGS) $(LIBS) 

plaquette_check: $(OBJECTS)
	$(CC) -o $@ $@.c $^ $(CFLAGS) $(WARNINGS) $(LIBS)

-include $(DEPENDS)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.c Makefile
	$(CC) -c $< -o $@  $(CFLAGS) $(WARNINGS) -MMD -MP

clean:
	rm -f $(BINARIES) $(BUILD_DIR)/*.[od] *~ */*~ 