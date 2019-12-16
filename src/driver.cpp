/*
 * driver.cpp
 */

#include "physics.h"
#include "voxels.h"

#define TEST_POP_SIZE 1

int main(int argc, char** argv) {


    Voxel_space* population[TEST_POP_SIZE];
    for (int i=0; i<TEST_POP_SIZE; i++) {
        INIT_VOXEL_SPACE(indiv);
        population[i] = indiv;
    }

    simulate_population_cpu(population, TEST_POP_SIZE, DEFAULT_START_HEIGHT);

    for (int i=0; i<TEST_POP_SIZE; i++) {
        printf("indiv %d distance travelled: %f, fitness: %f\n",
               i,
               population[i]->simulated_dist, 
               population[i]->fitness);
        export_to_gl(population[i], DEFAULT_START_HEIGHT);
    }

    for (int i=0; i<TEST_POP_SIZE; i++) {
        delete_voxel_space(population[i]);
    }
}

