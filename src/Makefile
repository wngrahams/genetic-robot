CC  = gcc
CXX = g++

INCLUDES = 

CFLAGS   = -O0 -g -Wall $(INCLUDES)
CXXFLAGS = -O0 -g -Wall $(INCLUDES)

LDFLAGS = -g 
LDLIBS  = -lm 

executables = driver
objects = driver.o ga-utils.o voxels.o physics.o

.PHONY: default
default: $(executables)

$(executables): ga-utils.o voxels.o physics.o

$(objects): ga-utils.h voxels.h physics.h

.PHONY: clean
clean:
	rm -f *~ a.out gl.log core $(objects) $(executables)

.PHONY: all
all: clean default
