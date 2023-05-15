INCLUDE_DIR := include
SOURCE_DIR := src
BUILD_DIR := build
BINS_DIR := bins

RM = rm -f
CC := gcc
CFLAGS := -I$(INCLUDE_DIR) -std=c99 -O4 -march=native -mtune=native -fopenmp -w -DNUM_THREADS=`nproc` -DAVX2

WARNINGS:= -Wall -Wextra -Werror

LIBS := -lm


SOURCES := $(wildcard $(SOURCE_DIR)/*.c)

_OBJECTS = $(patsubst $(SOURCE_DIR)/%.c, %.o, $(SOURCES))
OBJECTS = $(patsubst %, $(BUILD_DIR)/%, $(_OBJECTS))

_DEPENDS = $(patsubst $(SOURCE_DIR)/%.c, %.d, $(SOURCES))
DEPENDS = $(patsubst %, $(BUILD_DIR)/%, $(_DEPENDS))

BINARIES = gfix configGen invertDirac

.PHONY: all clean

all: $(BINARIES)


gfix: $(OBJECTS) gfix.c
	$(CC) -o $(BINS_DIR)/$@ $^ $(CFLAGS) $(WARNINGS) $(LIBS) 

configGen: $(OBJECTS) configGen.c
	$(CC) -o $(BINS_DIR)/$@ $^ $(CFLAGS) $(WARNINGS) $(LIBS) 


invertDirac: $(OBJECTS) invertDirac.c
	$(CC) -o $(BINS_DIR)/$@ $^ $(CFLAGS) $(WARNINGS) $(LIBS) 

-include $(DEPENDS)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.c Makefile
	$(CC) -c $< -o $@  $(CFLAGS) $(WARNINGS) -MMD -MP

clean:
	$(RM) $(BINARIES) $(OBJECTS) $(DEPENDS)