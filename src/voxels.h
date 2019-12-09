/*
 * voxels.h
 *
 */

#ifndef __VOXELS_H__
#define __VOXELS_H__

#include "ga-utils.h"

#define DEPTH_3X3   1
#define DEPTH_5X5   2
#define DEPTH_7X7   3
#define DEPTH_9X9   4
#define DEPTH_11X11 5

#define VOX_SPACE_MAX_DEPTH DEPTH_5X5  // Must be at least (3x3)

#define VOX_ROOT   0
#define VOX_MIDDLE 1
#define VOX_EDGE   2
#define VOX_CORNER 3

#define NUM_R_3x3 1
#define NUM_M_3X3 6
#define NUM_E_3X3 12
#define NUM_C_3X3 8


/*
 * struct Voxel: contains type indicator and position in voxel space
 */
typedef struct Voxel {

    int type;
    int pos_x;
    int pos_y;
    int pos_z;
    int exists;
    
} __attribute__ ((aligned)) Voxel;

/*
 * struct Voxel_space: contains array of voxels and the overall fitness
 */
typedef struct Voxel_space {

    Voxel* tree;
    float fitness;

} __attribute__ ((aligned)) Voxel_space;

/*
 * given a depth of a tree, returns the total number of voxels in that tree
 */
static inline int total_voxels_at_depth(int depth) {
    return ipow(2*depth + 1, 3);
}

static inline int num_root_at_depth(int depth) { return 1; }
int num_middle_at_depth(int);
int num_edge_at_depth(int);
static inline int num_corner_at_depth(int depth) { return 8*depth; }

void init_voxel_space(Voxel_space*);
void init_3x3(Voxel*);
void add_remaining_voxels(Voxel*);

#endif

