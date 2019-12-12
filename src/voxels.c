/*
 * voxels.c
 *
 */


#include "voxels.h"

/*
 * returns the total number of 'middle' voxels in a tree with max depth d
 */
int num_middle_at_depth(const int d) {
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

/*
 * returns the total number of 'edge' voxels in a tree with max depth d
 */
int num_edge_at_depth(const int d) {
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

/*
 * Initializes the entire voxel space over VOX_SPACE_MAX_DEPTH.
 * Allocates memory for the Voxel_space struct and the Voxel* array it contains.
 * Initializes all voxels' exists = 0 except for the root voxel.
 * Initializes all materials to UNKNOWN.
 */
void init_voxel_space(Voxel_space* vs) {

    //vs = malloc(sizeof(Voxel_space));
    //CHECK_MALLOC_ERR(vs);

    vs->tree = malloc(sizeof(Voxel)*total_voxels_at_depth(VOX_SPACE_MAX_DEPTH));
    CHECK_MALLOC_ERR(vs->tree);

    init_3x3(vs->tree);
    for (int i=0; i<total_voxels_at_depth(VOX_SPACE_MAX_DEPTH-1); i++) {
        init_children_of_index(vs->tree, i);
    }

    vs->num_voxels = 1;
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
    tree[0].type = ROOT;
    tree[0].position[0] = 0; 
    tree[0].position[1] = 0; 
    tree[0].position[2] = 0;
    tree[0].exists = 1;
    tree[0].material = UNKNOWN;

    // 'middle' voxels change the position in one dimension:
    start = NUM_R_3X3;
    stop = NUM_R_3X3 + NUM_M_3X3;
    for (int i=0; i<(stop-start); i++) {
        tree[i+start].type = MIDDLE;
        tree[i+start].position[0] = m_positions[3*i+0];
        tree[i+start].position[1] = m_positions[3*i+1];
        tree[i+start].position[2] = m_positions[3*i+2];
        tree[i+start].exists = 0;
        tree[i+start].material = UNKNOWN;
    }
    
    // 'edge' voxels change the position in two dimensions:
    start = stop;
    stop += NUM_E_3X3;
    for (int i=0; i<(stop-start); i++) {
        tree[i+start].type = EDGE;
        tree[i+start].position[0] = e_positions[3*i+0];
        tree[i+start].position[1] = e_positions[3*i+1];
        tree[i+start].position[2] = e_positions[3*i+2];
        tree[i+start].exists = 0;
        tree[i+start].material = UNKNOWN;
    }

    // 'corner' voxels change the position in three dimensions
    start = stop;
    stop += NUM_C_3X3;
    for (int i=0; i<(stop-start); i++) {
        tree[i+start].type = CORNER;
        tree[i+start].position[0] = c_positions[3*i+0];
        tree[i+start].position[1] = c_positions[3*i+1];
        tree[i+start].position[2] = c_positions[3*i+2];
        tree[i+start].exists = 0;
        tree[i+start].material = UNKNOWN;
    }

}

/*
 * initialize the childern of the voxel at index parent_idx
 */
void init_children_of_index(Voxel* tree, const int parent_idx) {
    
    switch(tree[parent_idx].type) {

        case CORNER:
            init_c_children(tree, parent_idx);
            break;

        case EDGE:
            init_e_children(tree, parent_idx);
            break;

        case MIDDLE:
            init_m_child(tree, parent_idx);
            break;

        case ROOT:
            break;
    }

}

/*
 * initializes the children of a 'corner' voxel at index parent_idx
 */
void init_c_children(Voxel* tree, const int parent_idx) {

    int parent_depth = get_depth_from_index(parent_idx);

    // 1 c child:
    int c_idx = total_voxels_at_depth(parent_depth+1)
                - total_voxels_at_depth(parent_depth)
                + parent_idx;

    tree[c_idx].type = CORNER;

    for (int i=0; i<3; i++) {
        tree[c_idx].position[i] = get_sign(tree[parent_idx].position[i])
                                  + tree[parent_idx].position[i];
    }

    tree[c_idx].exists = 0;
    tree[c_idx].material = UNKNOWN;

    // 3 m children:
    int* m_children = malloc(sizeof(int) * 3);
    CHECK_MALLOC_ERR(m_children);
    for (int i=0; i<3; i++) {
        m_children[i] = total_voxels_at_depth(parent_depth)
                      + num_middle_at_depth(parent_depth + 1) 
                      - num_middle_at_depth(parent_depth)
                      - 3*(total_voxels_at_depth(parent_depth) - parent_idx)
                      + i;
        tree[m_children[i]].type = MIDDLE;
        tree[m_children[i]].exists = 0;
        tree[m_children[i]].material = UNKNOWN;
    }
    init_m_positions(tree, parent_idx, m_children, 3);
    free(m_children);

    // 3 e children:
    int* e_children = malloc(sizeof(int) * 3);
    CHECK_MALLOC_ERR(e_children);
    for (int i=0; i<3; i++) {
        e_children[i] = total_voxels_at_depth(parent_depth)
                        + num_middle_at_depth(parent_depth + 1)
                        - num_middle_at_depth(parent_depth)
                        + num_edge_at_depth(parent_depth + 1)
                        - num_edge_at_depth(parent_depth)
                        - 3*(total_voxels_at_depth(parent_depth) - parent_idx)
                        + i;

        tree[e_children[i]].type = EDGE;
        tree[e_children[i]].exists = 0;
        tree[e_children[i]].material = UNKNOWN;
    }
    init_e_positions(tree, parent_idx, e_children, 3);
    free(e_children);
}

/*
 * initializes the children of an 'edge' voxel at parent_idx
 */
void init_e_children(Voxel* tree, const int parent_idx) {
    
    int parent_depth = get_depth_from_index(parent_idx);

    // 2 m children:
    int* m_children = malloc(sizeof(int) * 2);
    CHECK_MALLOC_ERR(m_children);
    for (int i=0; i<2; i++) {
        m_children[i] = total_voxels_at_depth(parent_depth)
                        + num_middle_at_depth(parent_depth)
                        - num_middle_at_depth(parent_depth - 1)
                        + 2*(parent_idx 
                            - (num_middle_at_depth(parent_depth)
                                - num_middle_at_depth(parent_depth - 1))
                            - total_voxels_at_depth(parent_depth - 1))
                        + i;
        tree[m_children[i]].type = MIDDLE;
        tree[m_children[i]].exists = 0;
        tree[m_children[i]].material = UNKNOWN;
    }
    init_m_positions(tree, parent_idx, m_children, 2);
    free(m_children);

    // 1 e child:
    int* e_children = malloc(sizeof(int));
    CHECK_MALLOC_ERR(e_children);
    for (int i=0; i<1; i++) {
        e_children[i] = total_voxels_at_depth(parent_depth)
                        + num_middle_at_depth(parent_depth + 1)
                        - num_middle_at_depth(parent_depth)
                        + parent_idx
                        - (num_middle_at_depth(parent_depth)
                            - num_middle_at_depth(parent_depth - 1))
                        - total_voxels_at_depth(parent_depth - 1);
        
        tree[e_children[i]].type = EDGE;
        tree[e_children[i]].exists = 0;
        tree[e_children[i]].material = UNKNOWN;
    }
    init_e_positions(tree, parent_idx, e_children, 1);
    free(e_children);
}

void init_m_child(Voxel* tree, const int parent_idx) {

    int parent_depth = get_depth_from_index(parent_idx);
    
    // 1 m child:
    int* m_child = malloc(sizeof(int));
    *m_child = total_voxels_at_depth(parent_depth)
               - total_voxels_at_depth(parent_depth - 1)
               + parent_idx;
    tree[*m_child].type = MIDDLE;
    tree[*m_child].exists = 0;
    tree[*m_child].material = UNKNOWN;

    init_m_positions(tree, parent_idx, m_child, 1);
    free(m_child);
}

/*
 * initializes the position member of the given indicies of 'middle' voxels
 * in the voxel space
 */
void init_m_positions(Voxel* tree, 
                      const int parent_idx, 
                      int* children_indices, 
                      const int num_children) {

    int* sorted_indices = malloc(sizeof(int) * 3);
    CHECK_MALLOC_ERR(sorted_indices);
    get_sorted_indices(tree[parent_idx].position, sorted_indices);

    for (int i=0; i<num_children; i++) {
        for (int j=0; j<3; j++) {
            tree[children_indices[i]].position[j] = tree[parent_idx].position[j];
        }

        tree[children_indices[i]].position[sorted_indices[i]]
            += get_sign(tree[parent_idx].position[sorted_indices[i]]);
    }

    free(sorted_indices);
}

/*
 * initializes the position member of the given indices of 'edge' voxels in the
 * voxel space
 */
void init_e_positions(Voxel* tree,
                      const int parent_idx,
                      int* children_indices,
                      const int num_children) {

    int* sorted_indices = malloc(sizeof(int) * 3);
    CHECK_MALLOC_ERR(sorted_indices);
    get_sorted_indices(tree[parent_idx].position, sorted_indices);

    for (int i=0; i<num_children; i++) {
        for (int j=0; j<3; j++) {
            tree[children_indices[i]].position[j] = tree[parent_idx].position[j];
        }

        tree[children_indices[i]].position[sorted_indices[i]]
            += get_sign(tree[parent_idx].position[sorted_indices[i]]);

        tree[children_indices[i]].position[sorted_indices[MOD(i+1,3)]]
            += get_sign(tree[parent_idx].position[sorted_indices[MOD(i+1,3)]]);
    }

    free(sorted_indices);
}

/*
 * returns the indeces of the position array in decreasing order of the 
 * absolute value of the position in each dimension
 */
void get_sorted_indices(int* position, int* sorted_indices) {

    int pos_abs[3] = { abs(position[0]), abs(position[1]), abs(position[2]) };

    if (pos_abs[0] == pos_abs[1] && pos_abs[0] == pos_abs[2]) {
        for (int i=0; i<3; i++) {
            sorted_indices[i] = i;
        }
    }
    else {
        int max_idx = 0;
        int min_idx = 0;
        int mid_idx = 0;

        for (int i=0; i<3; i++) {
            if (pos_abs[i] > pos_abs[max_idx])
                max_idx = i;
            if (pos_abs[i] < pos_abs[min_idx])
                min_idx = i;
        }

        mid_idx = 3 - (max_idx + min_idx);

        sorted_indices[0] = max_idx;
        sorted_indices[1] = mid_idx;
        sorted_indices[2] = min_idx;
    }
}

/*
 * free all memory associated with the voxel space
 */
void delete_voxel_space(Voxel_space* vs) {

    free(vs->tree);
    free(vs);
}

