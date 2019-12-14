/*
 * driver.c
 */

#include "physics.h"
#include "voxels.h"

#define TEST_POP_SIZE 4

int main(int argc, char** argv) {


    Voxel_space* population[TEST_POP_SIZE];
    for (int i=0; i<TEST_POP_SIZE; i++) {
        INIT_VOXEL_SPACE(indiv);
        population[i] = indiv;
    }

    simulate_population_cpu(population, TEST_POP_SIZE);

    for (int i=0; i<TEST_POP_SIZE; i++) {
        delete_voxel_space(population[i]);
    }
}

