#ifndef _GA_H
#define _GA_H

#include "voxels.h"

// GA macros

#define POP_SIZE 30 // should be even
#define NUM_OF_EVALS 240 // ideally a multiple of POP_SIZE
#define CHANCE_OF_MUT 0.5
#define NUM_OF_MUT 2

// robot macros
#define NUM_OF_CENTERS 5
#define NUM_OF_HOLES 9

#define NUM_OF_MATERIALS 4
#define NUM_OF_M_CHILD 1
#define NUM_OF_E_CHILD 3
#define NUM_OF_C_CHILD 7

// file macros
#define LEARNING_TXT "learningcurve_ga.txt" // for learning curve
#define DOT_TXT "dotplot_ga.txt"
#define DIVERSITY_TXT "diversity_ga.txt"

#define START_HEIGHT DEFAULT_START_HEIGHT

void ga_loop(int);
void initialize_random_robot(Voxel_space *);
void update_mats(Voxel_space *, const int, material_t);
void update_exists(Voxel_space* , const int, int);
void voxel_space_copy(Voxel_space *, Voxel_space *);
void crossover(Voxel_space *, Voxel_space *);
void crossover_exists(Voxel_space *, Voxel_space *, int);
void crossover_mats(Voxel_space *, Voxel_space *, int);
void mutation(Voxel_space *);
void selection(Voxel_space **, Voxel_space **, Voxel_space **);
void randomize_array(int *);
void copy_vs(Voxel_space *, Voxel_space *);
double calculate_diversity(Voxel_space **);
void random_loop(int);
void hc_loop(int);

#endif //_GA_H

