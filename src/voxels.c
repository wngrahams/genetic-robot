/*
 * voxels.c
 *
 */

#include "voxels.h"

int num_middle_at_depth(int d) {
    // base cases:
    if (0 == d)
        return 0;
    else if (1 == d)
        return NUM_M_3X3;

    // cached values:
    else if (2 == d)
        return 60;
    else if (3 == d)
        return 210;

    // recursive case:
    else {
        return 2*num_middle_at_depth(d-1) - num_middle_at_depth(d-2) 
                + 2*(num_edge_at_depth(d-1) - num_edge_at_depth(d-2)) + 24;
    }
}

int num_edge_at_depth(int d) {
    // base cases:
    if (0 == d)
        return 0;
    else if (1 == d)
        return NUM_E_3X3;

    // cached values:
    else if (2 == d)
        return 48;
    else if (3 == d)
        return 108;

    // recursive case:
    else 
        return 2*num_edge_at_depth(d-1) - num_edge_at_depth(d-2) + 24;
}

void init_voxel_space(Voxel_space* vs) {
    
    vs = malloc(sizeof(Voxel_space));
    CHECK_MALLOC_ERR(vs);

    vs->tree = malloc(sizeof(Voxel)*total_voxels_at_depth(VOX_SPACE_MAX_DEPTH));
    CHECK_MALLOC_ERR(vs->tree);

    init_3x3(vs->tree);
    add_remaining_voxels(vs->tree);
}

/*
 * initializes the inner 3x3 of the voxel space in such a way that the voxels 
 * are ordered by type: m->e->c and there is tight linkage between neighboring
 * children
 */
void init_3x3(Voxel* tree) {

    int start, stop;
    int m_positions[] = { 0,0,1, 1,0,0, 0,1,0, -1,0,0, 0,-1,0, 0,0,-1 };
    int e_positions[] = { -1,0,-1, 0,1,-1, 1,0,-1, 0,-1,-1, -1,-1,0, -1,1,0,
                          1,1,0, 1,-1,0, 0,-1,1, -1,0,1, 0,1,1, 1,0,1 };
    int c_positions[] = { 1,-1,1, -1,-1,1, -1,1,1, 1,1,1, 
                          1,-1,-1, -1,-1,-1, -1,1,-1, 1,1,-1 };

    // 'root' voxel at 0, 0, 0
    tree[0].type = VOX_ROOT;
    tree[0].pos_x = 0; tree[0].pos_y = 0; tree[0].pos_z = 0;
    tree[0].exists = 1;

    // 'middle' voxels change the position in one dimension:
    start = NUM_R_3X3;
    stop = NUM_R_3X3 + NUM_M_3X3;
    for (int i=start; i<stop; i++) {
        tree[i].type = VOX_MIDDLE;
        tree[i].pos_x = m_positions[3*i+0];
        tree[i].pos_y = m_positions[3*i+1];
        tree[i].pos_z = m_positions[3*i+2];
        tree[i].exists = 0;
    }
    
    // 'edge' voxels change the position in two dimensions:
    start = stop;
    stop += NUM_E_3X3;
    for (int i=start; i<stop; i++) {
        tree[i].type = VOX_EDGE;
        tree[i].pos_x = e_positions[3*i+0];
        tree[i].pos_y = e_positions[3*i+1];
        tree[i].pos_z = e_positions[3*i+2];
        tree[i].exists = 0;
    }

    // 'corner' voxels change the position in three dimensions
    start = stop;
    stop += NUM_C_3X3;
    for (int i=start; i<stop; i++) {
        tree[i].type = VOX_CORNER;
        tree[i].pos_x = c_positions[3*i+0];
        tree[i].pos_y = c_positions[3*i+1];
        tree[i].pos_z = c_positions[3*i+2];
        tree[1].exists = 0;
    }
}

/*
 * Fills in remainaing voxels with exists=0 for the total depth of the voxel 
 * space
 */
void add_remaining_voxels(Voxel* tree) {
    
}

