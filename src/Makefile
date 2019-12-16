CC  = gcc
CXX = g++

INCLUDES =

CFLAGS   = -O0 -g -Wall $(INCLUDES)
CXXFLAGS = -O0 -g -Wall -std=c++11 $(INCLUDES)

LDFLAGS = -g  
LDLIBS  = -lm -lstdc++ 

executables = driver
objects = drive.o ga-utils.o voxels.o physics.o obj_parser.o

.PHONY: default
default: $(executables)

$(executables): ga-utils.o voxels.o physics.o obj_parser.o

$(objects): ga-utils.h voxels.h physics.h obj_parser.hpp

.PHONY: clean
clean:
	rm -f *~ a.out gl.log core $(objects) $(executables)

.PHONY: all
all: clean default

