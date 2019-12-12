/*
 * driver.c
 */

#include "voxels.h"

int main(int argc, char** argv) {

    INIT_VOXEL_SPACE(indiv);

    /*
    for (int i=0; i<total_voxels_at_depth(VOX_SPACE_MAX_DEPTH); i++) {
        get_depth_from_index(i);
    }*/

    
    for (int i=0; i<total_voxels_at_depth(VOX_SPACE_MAX_DEPTH); i++) {
        printf("%d: %d (%d, %d, %d)\n", i, 
                                        indiv->tree[i].type,
                                        indiv->tree[i].position[0],
                                        indiv->tree[i].position[1],
                                        indiv->tree[i].position[2] );
        

        printf("%d: ", i);
        printf("%d ", indiv->tree[i].type);
        printf("(%d, ", indiv->tree[i].position[0]);
        printf("%d, ", indiv->tree[i].position[1]);
        printf("%d)\n", indiv->tree[i].position[2]);
    }
    delete_voxel_space(indiv);
}

