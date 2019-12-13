/*
 * voxels.h
 *
 */

#ifndef __VOXELS_H__
#define __VOXELS_H__

#include <math.h>
#include <stdlib.h>

#include "ga-utils.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DEPTH_3X3   1
#define DEPTH_5X5   2
#define DEPTH_7X7   3
#define DEPTH_9X9   4
#define DEPTH_11X11 5

#define VOX_SPACE_MAX_DEPTH DEPTH_11X11  // Must be at least (3x3)

#define NUM_R_3X3 1
#define NUM_M_3X3 6
#define NUM_E_3X3 12
#define NUM_C_3X3 8

#define INIT_VOXEL_SPACE(name) \
    Voxel_space* name = (Voxel_space*)malloc(sizeof(Voxel_space)); \
    CHECK_MALLOC_ERR((name)); \
    init_voxel_space((name)); 


typedef enum voxel_type {

    ROOT,
    MIDDLE,
    EDGE,
    CORNER

} __attribute__ ((packed)) voxel_type;

typedef enum material { 

    UNKNOWN, 
    HARD, 
    SOFT, 
    EXPAND, 
    CONTRACT 

} __attribute__ ((packed)) material_t;

/*
 * struct Voxel: contains type indicator and position in voxel space
 */
typedef struct Voxel {

    voxel_type type;
    int position[3];
    int exists;
    material_t material;
    
} Voxel;

/*
 * struct Voxel_space: contains array of voxels and the overall fitness
 */
typedef struct Voxel_space {

    Voxel* tree;
    int num_voxels;
    float fitness;

} Voxel_space;

/*
 * given a depth of a tree, returns the total number of voxels in that tree
 */
static inline int total_voxels_at_depth(const int depth) {
    return ipow(2*depth + 1, 3);
}

/*
 * given the index of a voxel, returns its depth in the tree
 */
static inline int get_depth_from_index(const int idx) {

    // this sucks, but without adding 0.000001, casting 5.0 to and int returns 
    // 4 for some reason
    int retval = (idx == 0) ? 0 : ceil(((int)((pow(idx, 1/3.)+0.000001))+0.0)/2.);
    return retval; 
}

/*
 * returns the total number of 'root' voxels in a tree with max depth d
 */
static inline int num_root_at_depth(const int d) { return 1; }

int num_middle_at_depth(const int);
int num_edge_at_depth(const int);

/*
 * returns the total number of 'corner' voxels in a tree with max depth d
 */
static inline int num_corner_at_depth(const int d) { return 8*d; }

void init_voxel_space(Voxel_space*);
void init_3x3(Voxel*);
void init_children_of_index(Voxel*, const int);
void init_c_children(Voxel*, const int);
void init_e_children(Voxel*, const int);
void init_m_child   (Voxel*, const int);
void init_m_positions(Voxel*, const int, int*, const int);
void init_e_positions(Voxel*, const int, int*, const int);

void get_sorted_indices(int*, int*);

void delete_voxel_space(Voxel_space*);

int get_child_index_of_m(const int);
int get_child_index_of_e(const int, const voxel_type, const int);
int get_child_index_of_c(const int, const voxel_type, const int);

#ifdef __cplusplus
}
#endif

#endif

