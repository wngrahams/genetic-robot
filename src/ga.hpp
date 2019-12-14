#ifndef _GA_H
#define _GA_H

#include "voxels.h"

// GA macros
#define POP_SIZE 4 // should be even
#define NUM_OF_EVALS 16 // ideally a multiple of POP_SIZE
#define CHANCE_OF_MUT 0.3
#define NUM_OF_MUT 1

// robot macros
#define NUM_OF_CENTERS 1
#define NUM_OF_HOLES 1

#define NUM_OF_MATERIALS 4
#define NUM_OF_M_CHILD 1
#define NUM_OF_E_CHILD 3
#define NUM_OF_C_CHILD 7

// file macros
#define LEARNING_TXT "learningcurve_ga.txt" // for learning curve

void ga_loop(int);
void initialize_random_robot(Voxel_space *);
void update_mats(Voxel_space *, const int, material_t);
void update_exists(Voxel_space* , const int, int);
void voxel_space_copy(Voxel_space *, Voxel_space *);
void crossover(Voxel_space *, Voxel_space *);
void crossover_exists(Voxel_space *, Voxel_space *, int p);
void mutation(Voxel_space *);
void selection(Voxel_space **, Voxel_space **, Voxel_space **);
void randomize_array(int *);

#endif //_GA_H