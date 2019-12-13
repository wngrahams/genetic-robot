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
#define L0_SIDE   0.1f
#define L0_FACE   0.141421f
#define L0_MIDDLE 0.173205f

#define K_HARD   10000.0f
#define K_SOFT   500.0f
#define K_MUSCLE 2000.0f

// These values must be multiplied by l0 before assigning to the b parameter!!!
#define B_STATIC 0.0f
#define B_MUSCLE 0.5f

static const float material_to_k_map[4]={ K_HARD,   K_SOFT,   K_MUSCLE, K_MUSCLE };
static const float material_to_b_map[4]={ B_STATIC, B_STATIC, B_MUSCLE, B_MUSCLE };
static const float material_to_c_map[4]={ 0.0f,     0.0f,     0.0f,     (F_PI/2.)};

static const float length_map[4]={ 0, L0_SIDE, L0_FACE, L0_MIDDLE };

typedef struct Mass {
    
    float m;       // mass in kilograms
    float pos[3];  // position in meters
    float vel[3];  // velocity in meters per second
    float acc[3];  // acceleration in meters per seconds squared

    int material;

} __attribute__ ((aligned(128))) Mass;

typedef struct Spring {

    unsigned int m1;
    unsigned int m2;

    float k;
    float l0;
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
int get_total_possible_masses(const int);
int get_total_possible_springs(const int); 
void dfs_init_masses(Voxel_space*, 
                     const int, 
                     Mass**, 
                     const int, 
                     int*);
void init_springs(Spring**, Mass**, int*, const int, const int);

static inline float get_b_from_mat(const int mat, const int l0) {
    return l0*mat;
}

#ifdef __cplusplus
}
#endif

#endif

