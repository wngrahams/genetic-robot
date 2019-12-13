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
    int max_robot_size = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);
    // create initial parent population
    Voxel_space parent[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++) {
        INIT_VOXEL_SPACE(temp);
        parent[i] = *temp;
    }

    // randomize initial population
    for (int i = 0; i < POP_SIZE; i++) {
        printf("%d: %d\n", i, parent[0].tree[i].type);
        initialize_random_robot(&parent[i]);
    }

    // genetic algorithm loop
    for (int eval = 0; eval < NUM_OF_EVALS; eval++) {
        // crossover

        // mutation

        // selection

    }

//    for (int i = 0; i < POP_SIZE; i++) {
//        delete_voxel_space(&parent[i]);
//    }
}

/*
 * initialize a robot with random material centers and morphology
 */
void initialize_random_robot(Voxel_space *individual) {

    // number of nodes in tree in Voxel_space
    int max_robot_size = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);
    std::cout << "max vox: " << max_robot_size << "\n\n\n\n";

    for (int i = 0; i < max_robot_size; i++) {
        individual->tree[i].exists = 1;
        individual->tree[i].material = UNKNOWN;
    }

    for (int i = 0; i < NUM_OF_CENTERS; i++) {
        // pick random index and random material
        int center = rand() % max_robot_size;
        material_t mat = static_cast<material_t>(rand() % (NUM_OF_MATERIALS - 1) + 1);
        // update tree accordingly
        update_mats(individual, center, mat);
    }

    for (int i = 0; i < NUM_OF_HOLES; i++) {
        // pick random index (not original cube) and set all children to null
        int center = rand() % (max_robot_size - 1) + 1;
        update_exists(individual, center, 0);
    }
}

/*
 * update voxel and all of its children with specified type
 */
void update_mats(Voxel_space* vs, const int idx, material_t mat) {

    std::cout << idx << " type: " << vs->tree[idx].type << "\n";
    if (idx >= total_voxels_at_depth(VOX_SPACE_MAX_DEPTH)) { std::cout << "HEY\n"; return; };

    for (int i = 0; i < total_voxels_at_depth(VOX_SPACE_MAX_DEPTH); i++) {
        std::cout << "typE: " << vs->tree[i].type << "\n";
    }

    // for convenience
    voxel_type current_type = vs->tree[idx].type;
    // set material
    vs->tree[idx].material = mat;

    // recursive descent through children
    if (ROOT == current_type) {
        for (int i=1; i<26; i++) {
            std::cout << i << ": ROOT HIT\n";
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

    if (idx >= total_voxels_at_depth(VOX_SPACE_MAX_DEPTH)) { std::cout << "EXISTS DONE\n"; return; };

    // for convenience
    voxel_type current_type = vs->tree[idx].type;
    // set material
    vs->tree[idx].exists = exists;

    // recursive descent through children
    if (ROOT == current_type) {
        for (int i=1; i<26; i++) {
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