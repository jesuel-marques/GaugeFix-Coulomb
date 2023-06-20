CC := gcc

BINARIES := gfix configGen invertDirac correlator fourier_transform_propagator averagePolyakov

DRIVER_DIR := driver
INCLUDE_DIR := include
SOURCE_DIR := src
BUILD_DIR := build
BINARIES_DIR := bin

CFLAGS := -I$(INCLUDE_DIR) -std=c99 -O3 -march=native -mtune=native -fopenmp -DNUM_THREADS=`nproc` -DAVX2 -flto 
DEPEND_FLAGS = -MMD -MP

WARNINGS := -Wall -Wextra -Wpedantic
WARNINGS_EXTRA := -Werror

LIBS := -lm

SOURCES = $(wildcard $(SOURCE_DIR)/*.c)
DRIVER_SOURCES = $(wildcard $(DRIVER_DIR)/*.c)

_OBJECTS = $(patsubst $(SOURCE_DIR)/%.c, %.o, $(SOURCES))
OBJECTS = $(patsubst %, $(BUILD_DIR)/%, $(_OBJECTS))

_DEPENDS = $(patsubst $(SOURCE_DIR)/%.c, %.d, $(SOURCES))
DEPENDS = $(patsubst %, $(BUILD_DIR)/%, $(_DEPENDS))

_DRIVER_OBJECTS = $(patsubst $(DRIVER_DIR)/%.c, %.o, $(DRIVER_SOURCES))
DRIVER_OBJECTS = $(patsubst %, $(BUILD_DIR)/%, $(_DRIVER_OBJECTS))

_DRIVER_DEPENDS = $(patsubst $(DRIVER_DIR)/%.c, %.d, $(DRIVER_SOURCES))
DRIVER_DEPENDS = $(patsubst %, $(BUILD_DIR)/%, $(_DRIVER_DEPENDS))

RM := rm -f

.SUFFIXES: .c .o
.PHONY: all clean

all: $(BINARIES)

gfix: $(BINARIES_DIR)/gfix
configGen: $(BINARIES_DIR)/configGen
invertDirac: $(BINARIES_DIR)/invertDirac
correlator: $(BINARIES_DIR)/correlator
fourier_transform_propagator: $(BINARIES_DIR)/fourier_transform_propagator
averagePolyakov: $(BINARIES_DIR)/averagePolyakov

$(BINARIES_DIR)/gfix: $(BUILD_DIR)/gfix.o $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(WARNINGS) $(LIBS)

$(BINARIES_DIR)/configGen: $(BUILD_DIR)/configGen.o $(OBJECTS) 
	$(CC) -o $@ $^ $(CFLAGS) $(WARNINGS) $(LIBS) 

$(BINARIES_DIR)/invertDirac: $(BUILD_DIR)/invertDirac.o $(OBJECTS) 
	$(CC) -o $@ $^ $(CFLAGS) $(WARNINGS) $(LIBS) 

$(BINARIES_DIR)/correlator: $(BUILD_DIR)/correlator.o $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(WARNINGS) $(LIBS) 

$(BINARIES_DIR)/fourier_transform_propagator: $(BUILD_DIR)/fourier_transform_propagator.o $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(WARNINGS) $(LIBS) 

$(BINARIES_DIR)/averagePolyakov: $(BUILD_DIR)/averagePolyakov.o $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(WARNINGS) $(LIBS) 	

-include $(DEPENDS) $(DRIVER_DEPENDS)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.c Makefile
	@if [ ! -d build ]; then mkdir build; fi
	$(CC) -c $< -o $@  $(CFLAGS) $(WARNINGS) $(DEPEND_FLAGS)

$(BUILD_DIR)/%.o: $(DRIVER_DIR)/%.c Makefile
	@if [ ! -d build ]; then mkdir build; fi
	$(CC) -c $< -o $@  $(CFLAGS) $(WARNINGS) $(DEPEND_FLAGS)

clean:
	$(RM) $(BINARIES) $(OBJECTS) $(DRIVER_OBJECTS) $(DEPENDS) $(DRIVER_DEPENDS)

rebuild: clean all