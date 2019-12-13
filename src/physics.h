/*
 * physics.h
 *
 */

#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include "ga-utils.h"
#include "voxels.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_MASSES_PER_VOXEL 8
#define MAX_SPRINGS_PER_VOXEL 28

#define HEIGHT_OFFSET VOX_SPACE_MAX_DEPTH

// cube variables
#define MASS_M 0.1f
#define L0_SIDE 0.1f

typedef struct Mass {
    
    float m;       // mass in kilograms
    float pos[3];  // position in meters
    float vel[3];  // velocity in meters per second
    float acc[3];  // acceleration in meters per seconds squared

} __attribute__ ((aligned(128))) Mass;

typedef struct Spring {

    float k;
    float l0;
    unsigned int m1;
    unsigned int m2;
    float a;
    float b;
    float c;

} __attribute__ ((aligned(64))) Spring;

void init_masses_and_springs_from_voxel_space(Mass**, 
                                              const int,
                                              Spring**,
                                              const int,
                                              int*,
                                              int*,
                                              Voxel_space*);
int get_total_possible_masses();
int get_total_possible_springs(); 
void dfs_init_masses(Voxel_space*, 
                     const int, 
                     Mass**, 
                     const int, 
                     int*);

#ifdef __cplusplus
}
#endif

#endif

