#include <iostream>
#include <fstream>
#include <random>
#include "ga.hpp"
#include "voxels.h"
#include "time.h"
#include "physics.h"

int main(int argc, char** argv) {

    srand(time(0));

//    INIT_VOXEL_SPACE(indiv);
//
//    for (int i = 0; i < total_voxels_at_depth(VOX_SPACE_MAX_DEPTH); i++) {
//        printf("%d: ", i);
//        printf("%d ", indiv->tree[i].type);
//        printf("(%d, ", indiv->tree[i].pos[0]);
//        printf("%d, ", indiv->tree[i].pos[1]);
//        printf("%d)\n", indiv->tree[i].pos[2]);
//    }
//    delete_voxel_space(indiv);

//    ga_loop(0);
    random_loop(0);
//    hc_loop(0);
}

/*
 * genetic algorithm loop
 */
void ga_loop(int thread_num) {

    // begin timer
    clock_t begin = clock();

    // initialize files
    std::ofstream learning_file;
    learning_file.open(std::to_string(thread_num) + LEARNING_TXT);
    std::ofstream dot_file;
    dot_file.open(std::to_string(thread_num) + DOT_TXT);
    std::ofstream diversity_file;
    diversity_file.open(std::to_string(thread_num) + DIVERSITY_TXT);

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
                   parent[i]->tree[j].pos[0],
                   parent[i]->tree[j].pos[1],
                   parent[i]->tree[j].pos[2],
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
                   parent[i]->tree[j].pos[0],
                   parent[i]->tree[j].pos[1],
                   parent[i]->tree[j].pos[2],
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
//            for (int j = 0; j < max_voxels; j++) {
//                memcpy(&child[i]->tree[j], &parent[i]->tree[j], sizeof(Voxel));
//            }
            copy_vs(child[i], parent[i]);
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
        for (int i = 0; i < POP_SIZE; i++) {
            simulate_population_cpu(child, POP_SIZE, START_HEIGHT);
        }

        // print fitnesses of population
        for (int i = 0; i < POP_SIZE; i++) {
            std::cout << 0 << ": " << parent[i]->fitness << "\n";
        }

        // initialize all population as copy of parent + child population
        Voxel_space* all[POP_SIZE * 2];
        for (int i = 0; i < POP_SIZE; i++) {
            INIT_VOXEL_SPACE(temp);
            all[i] = temp;
            all[i]->num_voxels = parent[i]->num_voxels;
            all[i]->fitness = parent[i]->fitness;
//            for (int j = 0; j < max_voxels; j++) {
//                memcpy(&all[i]->tree[j], &parent[i]->tree[j], sizeof(Voxel));
//            }
            copy_vs(all[i], parent[i]);
        }
        for (int i = POP_SIZE; i < POP_SIZE*2; i++) {
            INIT_VOXEL_SPACE(temp);
            all[i] = temp;
            all[i]->num_voxels = child[i - POP_SIZE]->num_voxels;
            all[i]->fitness = child[i - POP_SIZE]->fitness;
//            for (int j = 0; j < max_voxels; j++) {
//                memcpy(&all[i]->tree[j], &child[i - POP_SIZE]->tree[j], sizeof(Voxel));
//            }
            copy_vs(all[i], child[i - POP_SIZE]);
        }
        // selection
        selection(parent, child, all);

        // clean up child, all generation
        for (int i=0; i<POP_SIZE; i++) {
            delete_voxel_space(child[i]);
        }
        for (int i=0; i<POP_SIZE*2; i++) {
            delete_voxel_space(all[i]);
        }

        // find most fit individual
        int max_fit_index = 0;
        for (int i = 0; i < POP_SIZE; i++) {
            if (parent[i]->fitness > parent[max_fit_index]->fitness) {
                max_fit_index = i;
            }
        }

        // write learning curve to file
        for (int i = 0; i < POP_SIZE * 2; i++) {
            learning_file << parent[max_fit_index]->fitness << ",";
        }
        // write to dot plot file
        for (int i = 0; i < POP_SIZE; i++) {
            dot_file << eval << "," << parent[i]->fitness << "\n";
        }
        // calculate and write to diversity file
        diversity_file << calculate_diversity(parent) << ",";

        // print fitnesses of population
        if (eval % POP_SIZE == 0) {
            for (int i = 0; i < POP_SIZE; i++) {
                std::cout << eval << ": " << parent[i]->fitness << "\n";
            }
        }
    }
    
    // end timer
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    double iters_per_sec = NUM_OF_ITERATIONS / elapsed_secs;
    std::cout << "iter/sec: " << iters_per_sec << "\n";

    // find most fit individual
    int max_fit_index = 0;
    for (int i = 0; i < POP_SIZE; i++) {
        if (parent[i]->fitness > parent[max_fit_index]->fitness) {
            max_fit_index = i;
        }
    }
    std::cout << "most fit individual:\n";
    std::cout << "\tfitness: " << parent[max_fit_index]->fitness << std::endl;
    std::cout << "\tdist traveled: ";
    std::cout << parent[max_fit_index]->simulated_dist << std::endl;

    export_to_gl(parent[max_fit_index], START_HEIGHT);

    // clean up parent generation
    for (int i=0; i<POP_SIZE; i++) {
        delete_voxel_space(parent[i]);
    }

    // close file
    learning_file.close();
}

/*
 * initialize a robot with random material centers and morphology
 */
void initialize_random_robot(Voxel_space *individual) {
    // number of nodes in tree in Voxel_space
    int max_voxels = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);
    // random variable generation
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<> mat_init(0, NUM_MATERIALS-1);
    std::uniform_int_distribution<> mat_cent(0, max_voxels - 1);
    std::uniform_int_distribution<> exi_cent(1, max_voxels - 1);

    material_t start_mat = static_cast<material_t>(mat_init(mt));
    for (int i = 0; i < max_voxels; i++) {
        individual->tree[i].exists = 1;
        individual->tree[i].material = start_mat;
    }

    for (int i = 0; i < NUM_OF_CENTERS; i++) {
        // pick random index and a material
        int center = mat_cent(mt);
//        printf("center chosen: %d\n", center);
        material_t mat = static_cast<material_t>(mat_init(mt));
//        printf("material chosen: %d\n", mat);

        // update tree accordingly
        update_mats(individual, center, mat);
    }

    for (int i = 0; i < NUM_OF_HOLES; i++) {
        // pick random index (not root/central cube) and set all children to null
        int center = exi_cent(mt);
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
    int m = rand() % (max_voxels - 1) + 1;
    while (!A->tree[p].exists || !B->tree[p].exists) {
        p = rand() % (max_voxels - 1) + 1;
        m = rand() % (max_voxels - 1) + 1;
    }

    // switch morphology at point and all its children
    crossover_exists(A, B, p);
    crossover_mats(A, B, m);
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
 * recursively descend through two trees and switch exists at each idx
 */
void crossover_mats(Voxel_space *A, Voxel_space *B, int idx) {
    // return if past bounds
    if (idx >= total_voxels_at_depth(VOX_SPACE_MAX_DEPTH)) {
        return;
    }

    // for convenience
    voxel_type current_type = A->tree[idx].type;
    // switch morphology
    material_t temp = A->tree[idx].material;
    A->tree[idx].material = B->tree[idx].material;
    B->tree[idx].material = temp;

    // recursive descent through children
    if (ROOT == current_type) {
        for (int i=1; i<27; i++) {
            crossover_mats(A, B, i);
        }
    } else if (MIDDLE == current_type) {
        crossover_mats(A, B, get_child_index_of_m(idx));
    } else if (EDGE == current_type) {
        for (int i=0; i<2; i++) {
            crossover_mats(A, B, get_child_index_of_e(idx, MIDDLE, i));
        }
        crossover_mats(A, B, get_child_index_of_e(idx, EDGE, 0));
    } else if (CORNER == current_type) {
        for (int i=0; i<3; i++) {
            crossover_mats(A, B, get_child_index_of_c(idx, MIDDLE, i));
        }
        for (int i=0; i<3; i++) {
            crossover_mats(A, B, get_child_index_of_c(idx, EDGE, i));
        }
        crossover_mats(A, B, get_child_index_of_c(idx, CORNER, 0));
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
    for (int i = 0; i < NUM_OF_MUT; i++) {
        double mut = (rand()/(double)RAND_MAX);
        if (CHANCE_OF_MUT < mut) {
            // pick random spot for mutation
            int p = rand() % (max_voxels - 1) + 1;
            while (!A->tree[p].exists) {
                p = rand() % (max_voxels - 1) + 1;
            }
            // pick random exists or not and update robot
            material_t mat = (material_t) (rand() % NUM_OF_MATERIALS);
            update_mats(A, p, mat);
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
//    parent[0] = parent[best_parent_index];
//    parent[1] = child[best_child_index];
    copy_vs(parent[0], parent[best_parent_index]);
    copy_vs(parent[1], parent[best_child_index]);


    for (int i = 2; i < POP_SIZE; i++) {
        int m = compare(mt);
        int n = compare(mt);
        while (m == n) {
            n = compare(mt);
        }

        if (all[m]->fitness > all[n]->fitness) {
//            parent[i] = all[m];
            copy_vs(parent[i], all[m]);
        } else {
//            parent[i] = all[n];
            copy_vs(parent[i], all[n]);
        }
    }
}

/*
 * create random array
 */
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

/*
 * copy one Voxel_space to another
 */
void copy_vs(Voxel_space *child, Voxel_space *parent) {
    child->num_voxels = parent->num_voxels;
    child->fitness = parent->fitness;
    child->simulated_dist = parent->simulated_dist;
    for (int j = 0; j < total_voxels_at_depth(VOX_SPACE_MAX_DEPTH); j++) {
        child->tree[j].pos[0] = parent->tree[j].pos[0];
        child->tree[j].pos[1] = parent->tree[j].pos[1];
        child->tree[j].pos[2] = parent->tree[j].pos[2];
        child->tree[j].exists = parent->tree[j].exists;
        child->tree[j].material = parent->tree[j].material;
        child->tree[j].type = parent->tree[j].type;
    }
}

/*
 * calculate diversity based on avg position of empty cubes
 */
double calculate_diversity(Voxel_space **parent) {
    int max_voxels = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);
    int avg_pos = 0;
    double mse = 0.0;

    for (int i = 0; i < POP_SIZE; i++) {
        for (int j = 0; j < max_voxels; j++) {
            if (!parent[i]->tree[j].exists) {
                avg_pos += j;
            }
        }
    }
    avg_pos /= POP_SIZE * max_voxels;
    for (int i = 0; i < POP_SIZE; i++) {
        int sum = 0;
        for (int j = 0; j < max_voxels; j++) {
            if (!parent[i]->tree[j].exists) {
                sum += j;
            }
        }
        mse += pow(sum - avg_pos, 2);
    }
    mse /= POP_SIZE;
    return mse;
}

void random_loop(int thread_num) {

    // begin timer
    clock_t begin = clock();

    // initialize files
    std::ofstream learning_file;
    learning_file.open(std::to_string(thread_num) + LEARNING_TXT);
    std::ofstream dot_file;
    dot_file.open(std::to_string(thread_num) + DOT_TXT);
    std::ofstream diversity_file;
    diversity_file.open(std::to_string(thread_num) + DIVERSITY_TXT);

    // declare variables
    int max_voxels = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);

    // create initial parent population
    Voxel_space* parent[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++) {
        INIT_VOXEL_SPACE(temp);
        parent[i] = temp;
    }

    // randomize initial population
    for (int i = 0; i < POP_SIZE; i++) {
        //printf("%d: %d\n", i, parent[0].tree[i].type);
        initialize_random_robot(parent[i]);
    }

    // calculate fitness
    for (int i = 0; i < POP_SIZE; i++) {
        simulate_population_cpu(parent, POP_SIZE, START_HEIGHT);
    }

    // genetic algorithm loop
    for (int eval = 0; eval < NUM_OF_EVALS; eval+=POP_SIZE) {

        // create other parent population
        Voxel_space* other[POP_SIZE];
        for (int i = 0; i < POP_SIZE; i++) {
            INIT_VOXEL_SPACE(temp);
            other[i] = temp;
        }

        // randomize other population
        initialize_random_robot(other[0]);

        // calculate fitness
        for (int i = 0; i < POP_SIZE; i++) {
            simulate_population_cpu(parent, POP_SIZE, START_HEIGHT);
            simulate_population_cpu(other, POP_SIZE, START_HEIGHT);
        }

        // print fitnesses of population
        std::cout << 0 << ": " << parent[0]->fitness << "\n";

        // test
        if (parent[0]->fitness < other[0]->fitness) {
            copy_vs(parent[0], other[0]);
        }

        // write learning curve to file
        learning_file << parent[0]->fitness << ",";
        learning_file << parent[0]->fitness << ",";

        // print fitnesses of population
        if (eval % POP_SIZE == 0) {
            for (int i = 0; i < POP_SIZE; i++) {
                std::cout << eval << ": " << parent[i]->fitness << "\n";
            }
        }
    }

    // end timer
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    double iters_per_sec = NUM_OF_ITERATIONS / elapsed_secs;
    std::cout << "iter/sec: " << iters_per_sec << "\n";

    std::cout << "most fit individual:\n";
    std::cout << "\tfitness: " << parent[0]->fitness << std::endl;
    std::cout << "\tdist traveled: ";
    std::cout << parent[0]->simulated_dist << std::endl;

    export_to_gl(parent[0], START_HEIGHT);

    // clean up parent generation
    delete_voxel_space(parent[0]);

    // close file
    learning_file.close();
}

void hc_loop(int thread_num) {

    // begin timer
    clock_t begin = clock();

    // initialize files
    std::ofstream learning_file;
    learning_file.open(std::to_string(thread_num) + LEARNING_TXT);
    std::ofstream dot_file;
    dot_file.open(std::to_string(thread_num) + DOT_TXT);
    std::ofstream diversity_file;
    diversity_file.open(std::to_string(thread_num) + DIVERSITY_TXT);

    // declare variables
    int max_voxels = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);

    // create initial parent population
    Voxel_space* parent[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++) {
        INIT_VOXEL_SPACE(temp);
        parent[i] = temp;
    }
    // randomize initial population
    for (int i = 0; i < POP_SIZE; i++) {
        //printf("%d: %d\n", i, parent[0].tree[i].type);
        initialize_random_robot(parent[i]);
    }
    // calculate fitness
    for (int i = 0; i < POP_SIZE; i++) {
        simulate_population_cpu(parent, POP_SIZE, START_HEIGHT);
    }

    // hc loop
    for (int eval = 0; eval < NUM_OF_EVALS; eval+=POP_SIZE) {

        // create other parent population
        Voxel_space* other[POP_SIZE];
        for (int i = 0; i < POP_SIZE; i++) {
            INIT_VOXEL_SPACE(temp);
            other[i] = temp;
        }

        // randomize other population
        copy_vs(other[0], parent[0]);

        // mutation
        mutation(other[0]);

        // calculate fitness
        for (int i = 0; i < POP_SIZE; i++) {
            simulate_population_cpu(other, POP_SIZE, START_HEIGHT);
        }

        // print fitnesses of population
        std::cout << 0 << ": " << parent[0]->fitness << "\n";

        // test
        if (parent[0]->fitness < other[0]->fitness) {
            copy_vs(parent[0], other[0]);
        }

        // write learning curve to file
        learning_file << parent[0]->fitness << ",";
        learning_file << parent[0]->fitness << ",";

        // print fitnesses of population
        if (eval % POP_SIZE == 0) {
            for (int i = 0; i < POP_SIZE; i++) {
                std::cout << eval << ": " << parent[i]->fitness << "\n";
            }
        }
    }

    // end timer
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    double iters_per_sec = NUM_OF_ITERATIONS / elapsed_secs;
    std::cout << "iter/sec: " << iters_per_sec << "\n";

    std::cout << "most fit individual:\n";
    std::cout << "\tfitness: " << parent[0]->fitness << std::endl;
    std::cout << "\tdist traveled: ";
    std::cout << parent[0]->simulated_dist << std::endl;

    export_to_gl(parent[0], START_HEIGHT);

    // clean up parent generation
    delete_voxel_space(parent[0]);

    // close file
    learning_file.close();
}