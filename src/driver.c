/*
 * driver.c
 */

#include "physics.h"
#include "voxels.h"

int main(int argc, char** argv) {

    INIT_VOXEL_SPACE(indiv);

    for (int i=0; i<total_voxels_at_depth(VOX_SPACE_MAX_DEPTH); i++) {
        printf("%d: ", i);
        printf("%d ", indiv->tree[i].type);
        printf("(%d, ", indiv->tree[i].pos[0]);
        printf("%d, ", indiv->tree[i].pos[1]);
        printf("%d)\n", indiv->tree[i].pos[2]);
    }
    
    int max_masses = get_total_possible_masses(VOX_SPACE_MAX_DEPTH);
    Mass** possible_masses = malloc(sizeof(Mass*) * max_masses);

    int max_springs = get_total_possible_springs(VOX_SPACE_MAX_DEPTH);
    Spring** possible_springs = malloc(sizeof(Spring*) * max_springs);

    int num_actual_masses = 0;
    int num_actual_springs = 0;
    init_masses_and_springs_from_voxel_space(possible_masses,
                                             max_masses,
                                             possible_springs,
                                             max_springs,
                                             &num_actual_masses,
                                             &num_actual_springs,
                                             indiv);

    printf("num masses initialized: %d\n", num_actual_masses);
    printf("num springs initialized: %d\n", num_actual_springs);

    for (int i=0; i<max_masses; i++) {
        if (NULL != possible_masses[i]) {
            printf("mass %d with pos (%.2f, %.2f, %.2f)\n", i,
                   possible_masses[i]->pos[0],
                   possible_masses[i]->pos[1],
                   possible_masses[i]->pos[2]);

        }
    }

    for (int i=0; i<max_masses; i++) {
        if (NULL != possible_masses[i])
            free(possible_masses[i]);
    }
    free(possible_masses);

    for (int i=0; i<max_springs; i++) {
        if (NULL != possible_springs[i])
            free(possible_springs[i]);
    }
    free(possible_springs);



    delete_voxel_space(indiv);
}

