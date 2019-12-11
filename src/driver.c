/*
 * driver.c
 */

#include "voxels.h"

int main(int argc, char** argv) {

    INIT_VOXEL_SPACE(indiv);

    for (int i=0; i<indiv->num_voxels; i++) {
        printf("loop\n");
        printf("%d: %d (%d, %d, %d)\n", i, 
                                        indiv->tree[i].type,
                                        indiv->tree[i].position[0],
                                        indiv->tree[i].position[1],
                                        indiv->tree[i].position[2] );
    }

    delete_voxel_space(indiv);
}

