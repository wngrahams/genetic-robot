CC  = gcc-9
CXX = g++-9

INCLUDES = 

CFLAGS   = -O3 -g -Wall -Werror -fopenmp $(INCLUDES)
CXXFLAGS = -O3 -g -Wall -Werror -fopenmp -std=c++11 $(INCLUDES)

LDFLAGS = -g -fopenmp 
LDLIBS  = -lm -lstdc++ -pthread 

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
