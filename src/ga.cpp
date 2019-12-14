#include <iostream>
#include "ga.hpp"
#include "voxels.h"
#include "time.h"

int main(int argc, char** argv) {

    srand(time(0));

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

    ga_loop();
}

/*
 * genetic algorithm loop
 */
void ga_loop() {

    // declare variables
    int max_voxels = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);

    // create initial parent population
    Voxel_space* parent[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++) {
        INIT_VOXEL_SPACE(temp);
        parent[i] = temp;
    }

    for (int i=0; i<POP_SIZE; i++) {
        printf("individual %d:\n", i);
        for (int j=0; j<max_voxels; j++) {
            printf("\tvox[%d]: (%d, %d, %d), exists=%d, material=%d, type=%d\n",
                   j,
                   parent[i]->tree[j].position[0],
                   parent[i]->tree[j].position[1],
                   parent[i]->tree[j].position[2],
                   parent[i]->tree[j].exists,
                   parent[i]->tree[j].material,
                   parent[i]->tree[j].type);

        }
        printf("\n");
    }
    
    // randomize initial population
    for (int i = 0; i < POP_SIZE; i++) {
        //printf("%d: %d\n", i, parent[0].tree[i].type);
        initialize_random_robot(parent[i]);
    }

    for (int i=0; i<POP_SIZE; i++) {
        printf("parent %d:\n", i);
        for (int j=0; j<max_voxels; j++) {
            printf("\tvox[%d]: (%d, %d, %d), exists=%d, material=%d, type=%d\n",
                   j,
                   parent[i]->tree[j].position[0],
                   parent[i]->tree[j].position[1],
                   parent[i]->tree[j].position[2],
                   parent[i]->tree[j].exists,
                   parent[i]->tree[j].material,
                   parent[i]->tree[j].type);

        }
        printf("\n");
    }

    // genetic algorithm loop
    for (int eval = 0; eval < NUM_OF_EVALS; eval++) {
        // initialize children population as copy of parent population
        Voxel_space* child[POP_SIZE];
        for (int i = 0; i < POP_SIZE; i++) {
            INIT_VOXEL_SPACE(temp);
            child[i] = temp;
            child[i]->num_voxels = parent[i]->num_voxels;
            child[i]->fitness = parent[i]->fitness;
            for (int j = 0; j < max_voxels; j++) {
                memcpy(&child[i]->tree[j], &parent[i]->tree[j], sizeof(Voxel));
            }
        }

        for (int i=0; i<POP_SIZE; i++) {
            printf("child %d:\n", i);
            for (int j=0; j<max_voxels; j++) {
                printf("\tvox[%d]: (%d, %d, %d), exists=%d, material=%d, type=%d\n",
                       j,
                       child[i]->tree[j].position[0],
                       child[i]->tree[j].position[1],
                       child[i]->tree[j].position[2],
                       child[i]->tree[j].exists,
                       child[i]->tree[j].material,
                       child[i]->tree[j].type);

            }
            printf("\n");
        }

//        // crossover

//        // mutation

//        // selection

//        // clean up child generation
        for (int i=0; i<POP_SIZE; i++) {
            delete_voxel_space(child[i]);
        }

    }

//    for (int i = 0; i < POP_SIZE; i++) {
//        delete_voxel_space(&parent[i]);
//    }


    // clean up
    for (int i=0; i<POP_SIZE; i++) {
        delete_voxel_space(parent[i]);
    }
}

/*
 * initialize a robot with random material centers and morphology
 */
void initialize_random_robot(Voxel_space *individual) {

    // number of nodes in tree in Voxel_space
    int max_voxels = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);
    std::cout << "max vox: " << max_voxels << "\n";

    for (int i = 0; i < max_voxels; i++) {
        individual->tree[i].exists = 1;
        individual->tree[i].material = BONE;
    }

    for (int i = 0; i < NUM_OF_CENTERS; i++) {
        // pick random index and a material
        int center = rand() % max_voxels;
        printf("center chosen: %d\n", center);
        material_t mat = static_cast<material_t>(i);
        printf("material chosen: %d\n", mat);

        // update tree accordingly
        update_mats(individual, center, mat);
    }

    for (int i = 0; i < NUM_OF_HOLES; i++) {
        // pick random index (not root/central cube) and set all children to null
        int center = rand() % (max_voxels - 1) + 1;
        update_exists(individual, center, 0);
    }
}

/*
 * update voxel and all of its children with specified material
 */
void update_mats(Voxel_space* vs, const int idx, material_t mat) {
    // return if past bounds
    if (idx >= total_voxels_at_depth(VOX_SPACE_MAX_DEPTH)) {
        return;
    }

    // for convenience
    voxel_type current_type = vs->tree[idx].type;
    // set material
    vs->tree[idx].material = mat;

    // recursive descent through children
    if (ROOT == current_type) {
        for (int i=1; i<27; i++) {
            //std::cout << i << ": ROOT HIT\n";
            update_mats(vs, i, mat);
        }
    } else if (MIDDLE == current_type) {
        update_mats(vs, get_child_index_of_m(idx), mat);
    } else if (EDGE == current_type) {
        for (int i=0; i<2; i++) {
            update_mats(vs, get_child_index_of_e(idx, MIDDLE, i), mat);
        }
        update_mats(vs, get_child_index_of_e(idx, EDGE, 0), mat);
    } else if (CORNER == current_type) {
        for (int i=0; i<3; i++) {
            update_mats(vs, get_child_index_of_c(idx, MIDDLE, i), mat);
        }
        for (int i=0; i<3; i++) {
            update_mats(vs, get_child_index_of_c(idx, EDGE, i), mat);
        }
        update_mats(vs, get_child_index_of_c(idx, CORNER, 0), mat);
    }
}

/*
 * update voxel and all of its children with existence or non-existence
 */
void update_exists(Voxel_space* vs, const int idx, int exists) {
    // return if past bounds
    if (idx >= total_voxels_at_depth(VOX_SPACE_MAX_DEPTH)) {
        return;
    }

    // for convenience
    voxel_type current_type = vs->tree[idx].type;
    // set material
    vs->tree[idx].exists = exists;

    // recursive descent through children
    if (ROOT == current_type) {
        for (int i=1; i<27; i++) {
            update_exists(vs, i, exists);
        }
    } else if (MIDDLE == current_type) {
        update_exists(vs, get_child_index_of_m(idx), exists);
    } else if (EDGE == current_type) {
        for (int i=0; i<2; i++) {
            update_exists(vs, get_child_index_of_e(idx, MIDDLE, i), exists);
        }
        update_exists(vs, get_child_index_of_e(idx, EDGE, 0), exists);
    } else if (CORNER == current_type) {
        for (int i=0; i<3; i++) {
            update_exists(vs, get_child_index_of_c(idx, MIDDLE, i), exists);
        }
        for (int i=0; i<3; i++) {
            update_exists(vs, get_child_index_of_c(idx, EDGE, i), exists);
        }
        update_exists(vs, get_child_index_of_c(idx, CORNER, 0), exists);
    }
}

/*
 * copy one voxel_space array population to another
 */
void voxel_space_copy(Voxel_space *A, Voxel_space *B) {
    // implemented in ga_loop currently
}

/*
 * crossover
 */
void crossover(Voxel_space *vs) {

}

/*
 * mutation
 */
void mutation(Voxel_space *vs) {

}

/*
 * tournament selection
 */
void selection(Voxel_space *parent, Voxel_space *child) {

}
