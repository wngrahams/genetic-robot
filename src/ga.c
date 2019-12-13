//
// Created by Jarrett Ross on 12/12/19.
//

#include "ga.h"
#include "voxels.h"

int main(int argc, char** argv) {

//    INIT_VOXEL_SPACE(indiv);
//
//    for (int i = 0; i < total_voxels_at_depth(VOX_SPACE_MAX_DEPTH); i++) {
//        printf("%d: ", i);
//        printf("%d ", indiv->tree[i].type);
//        printf("(%d, ", indiv->tree[i].position[0]);
//        printf("%d, ", indiv->tree[i].position[1]);
//        printf("%d)\n", indiv->tree[i].position[2]);
//    }
//    delete_voxel_space(indiv);

    loop();
}

void loop() {

    Voxel_space parent[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++) {
        INIT_VOXEL_SPACE(temp);
        parent[i] = *temp;
    }

    for (int i = 0; i < total_voxels_at_depth(VOX_SPACE_MAX_DEPTH); i++) {
        printf("%d: %d\n", i, parent[0].tree[i].type);
    }

}

void initialize_random_robot(Voxel_space individual) {


    for (int i = 0; i <


}


