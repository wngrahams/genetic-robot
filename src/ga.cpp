#include <iostream>
#include "ga.hpp"
#include "voxels.h"
#include "time.h"
#include <random>

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
    for (int eval = 0; eval < NUM_OF_EVALS; eval+=POP_SIZE) {
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

        // generate random order for crossover
        int order[POP_SIZE];
        randomize_array(order);
        // crossover
        for (int i = 0; i < POP_SIZE; i+=2) {
            crossover(child[order[i]], child[order[i+1]]);
        }

        // mutation
        for (int i = 0; i < POP_SIZE; i++) {
            mutation(child[i]);
        }

        // calculate fitness
        // fitness this dicc in ur mouth

        // initialize all population as copy of parent + child population
        Voxel_space* all[POP_SIZE * 2];
        for (int i = 0; i < POP_SIZE; i++) {
            INIT_VOXEL_SPACE(temp);
            all[i] = temp;
            all[i]->num_voxels = parent[i]->num_voxels;
            all[i]->fitness = parent[i]->fitness;
            for (int j = 0; j < max_voxels; j++) {
                memcpy(&all[i]->tree[j], &parent[i]->tree[j], sizeof(Voxel));
            }
        }
        for (int i = POP_SIZE; i < POP_SIZE*2; i++) {
            INIT_VOXEL_SPACE(temp);
            all[i] = temp;
            all[i]->num_voxels = child[i]->num_voxels;
            all[i]->fitness = child[i]->fitness;
            for (int j = 0; j < max_voxels; j++) {
                memcpy(&all[i]->tree[j], &child[i]->tree[j], sizeof(Voxel));
            }
        }
        // selection
//        selection(parent, child, all);

        // clean up child, all generation
        for (int i=0; i<POP_SIZE; i++) {
            delete_voxel_space(child[i]);
        }
        for (int i=0; i<POP_SIZE*2; i++) {
            delete_voxel_space(all[i]);
        }
    }

    // clean up parent generation
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
void crossover(Voxel_space *A, Voxel_space *B) {
    // get random spot for crossover (not root / center)
    int max_voxels = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);
    int p = rand() % (max_voxels - 1) + 1;
    while (!A->tree[p].exists || !B->tree[p].exists) {
        p = rand() % (max_voxels - 1) + 1;
    }

    // switch morphology at point and all its children
    crossover_exists(A, B, p);
}

/*
 * recursively descend through two trees and switch exists at each idx
 */
void crossover_exists(Voxel_space *A, Voxel_space *B, int idx) {
    // return if past bounds
    if (idx >= total_voxels_at_depth(VOX_SPACE_MAX_DEPTH)) {
        return;
    }

    // for convenience
    voxel_type current_type = A->tree[idx].type;
    // switch morphology
    int temp = A->tree[idx].exists;
    A->tree[idx].exists = B->tree[idx].exists;
    B->tree[idx].exists = temp;

    // recursive descent through children
    if (ROOT == current_type) {
        for (int i=1; i<27; i++) {
            crossover_exists(A, B, i);
        }
    } else if (MIDDLE == current_type) {
        crossover_exists(A, B, get_child_index_of_m(idx));
    } else if (EDGE == current_type) {
        for (int i=0; i<2; i++) {
            crossover_exists(A, B, get_child_index_of_e(idx, MIDDLE, i));
        }
        crossover_exists(A, B, get_child_index_of_e(idx, EDGE, 0));
    } else if (CORNER == current_type) {
        for (int i=0; i<3; i++) {
            crossover_exists(A, B, get_child_index_of_c(idx, MIDDLE, i));
        }
        for (int i=0; i<3; i++) {
            crossover_exists(A, B, get_child_index_of_c(idx, EDGE, i));
        }
        crossover_exists(A, B, get_child_index_of_c(idx, CORNER, 0));
    }
}

/*
 * mutation
 */
void mutation(Voxel_space *A) {
    int max_voxels = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);
    for (int i = 0; i < NUM_OF_MUT; i++) {
        double mut = (rand()/(double)RAND_MAX);
        if (CHANCE_OF_MUT < mut) {
            // pick random spot for mutation
            int p = rand() % (max_voxels - 1) + 1;
            while (!A->tree[p].exists) {
                p = rand() % (max_voxels - 1) + 1;
            }
            // pick random exists or not and update robot
            int exists = rand() % 1;
            update_exists(A, p, exists);
        }
    }
}

/*
 * tournament selection
 */
void selection(Voxel_space **parent, Voxel_space **child, Voxel_space **all) {
    // random variable generation
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<> compare(0, POP_SIZE * 2 - 1);

    // take elite child;
    int best_parent_index = 0;
    int best_child_index = 0;
    for (int i = 0; i < POP_SIZE; i++) {
        if (parent[i]->fitness > parent[best_parent_index]->fitness) {
            best_parent_index = i;
        }
        if (child[i]->fitness > child[best_child_index]->fitness) {
            best_child_index = i;
        }
    }
    parent[0] = parent[best_parent_index];
    parent[1] = child[best_child_index];

    for (int i = 2; i < POP_SIZE; i++) {
        int m = compare(mt);
        int n = compare(mt);
        while (m == n) {
            n = compare(mt);
        }

        if (all[m]->fitness > all[n]->fitness) {
            parent[i] = all[m];
        } else {
            parent[i] = all[n];
        }
    }
}

void randomize_array(int *order) {
    // create vector with order of parents to be crossed over
    for (int i = 0; i < POP_SIZE; i++) {
        order[i] = i;
    }
    for (int i = POP_SIZE - 1; i > 0; i--) {
        int j = (int)(rand() % (i+1));
        int temp = order[i];
        order[i] = order[j];
        order[j] = temp;
    }
}
