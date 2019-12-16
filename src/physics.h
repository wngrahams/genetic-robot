/*
 * physics.h
 *
 */

#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include <assert.h>

#include "ga-utils.h"
#include "voxels.h"

#ifdef __cplusplus
extern "C" {
#endif

#define OBJ_FILE "vs.obj"

#define MAX_MASSES_PER_VOXEL 8
#define MAX_SPRINGS_PER_VOXEL 28

#define POS_OFFSET VOX_SPACE_MAX_DEPTH

// world variables
#define G                 9.80665f
#define OMEGA             (4.*F_PI)
#define K_GROUND          100000.0f
#define DT                0.0001f
#define V_DAMP_CONST      0.999f //0.999999
#define NUM_OF_ITERATIONS 15000
// 15000 = 3 cycles (for omega = 4pi and DT = 0.0001)
#define U_S 1.0f
#define U_K 0.8f

// cube variables
#define MASS_M    0.1f
#define L0_SIDE   0.1f
#define L0_FACE   0.141421f
#define L0_MIDDLE 0.173205f

#define DEFAULT_START_HEIGHT (L0_SIDE/100.0f)

#define K_HARD   10000.0f
#define K_SOFT   500.0f
#define K_MUSCLE 2000.0f

// These values must be multiplied by L0 before assigning to the b parameter!!!
#define B_STATIC 0.0f
#define B_MUSCLE 0.4f

static const float material_to_k_map[4]={ K_HARD,   K_SOFT,   K_MUSCLE, K_MUSCLE };
static const float material_to_b_map[4]={ B_STATIC, B_STATIC, B_MUSCLE, B_MUSCLE };
static const float material_to_c_map[4]={ 0.0f,     0.0f,     0.0f,     (F_PI/2.)};

static const float length_map[4]={ 0.0f, L0_SIDE, L0_FACE, L0_MIDDLE };

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
                                              Voxel_space*,
                                              const float);
int get_total_possible_masses(const int);
int get_total_possible_springs(const int); 
void dfs_init_masses(Voxel_space*, 
                     const int, 
                     Mass**, 
                     const int, 
                     int*,
                     const float);
void init_springs(Spring**, Mass**, int*, const int, const int);

static inline float get_b_from_mat(const int mat, const float l0) {
    assert(mat >=0 && mat < NUM_MATERIALS);
    return l0*material_to_b_map[mat];
}

void simulate_population_cpu(Voxel_space**, const int, const float); 

static float inline dist3d(const float x2, const float x1,
                           const float y2, const float y1,
                           const float z2, const float z1) {

    return sqrtf(powf(x2-x1, 2.) + powf(y2-y1, 2.) + powf(z2-z1, 2.));
}

void calculate_center_of_mass(Mass**, const int, float*);

void export_to_gl(Voxel_space*, const float);
void write_obj(Voxel_space*, Mass**, const int, const char*);
void simulate_gl(Voxel_space*, const float, char*);

#ifdef __cplusplus
}
#endif

#endif

