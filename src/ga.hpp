#ifndef _GA_H
#define _GA_H

#include "voxels.h"

// GA macros
#define POP_SIZE 4
#define NUM_OF_EVALS 1000

// robot macros
#define NUM_OF_CENTERS 1
#define NUM_OF_HOLES 1

#define NUM_OF_MATERIALS 4
#define NUM_OF_M_CHILD 1
#define NUM_OF_E_CHILD 3
#define NUM_OF_C_CHILD 7

// file macros
#define LEARNING_TXT "learningcurve_ga.txt" // for learning curve

void ga_loop();
void initialize_random_robot(Voxel_space *);
void update_mats(Voxel_space *, const int, material_t);
void update_exists(Voxel_space* , const int, int);

#endif //_GA_H
