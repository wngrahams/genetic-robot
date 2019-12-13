/*
 * physics.c
 *
 */

#include "physics.h"

void init_masses_and_springs_from_voxel_space(Mass** masses, 
                                              const int max_num_masses,
                                              Spring** springs,
                                              const int max_num_springs,
                                              int* mass_count,
                                              int* spring_count,
                                              Voxel_space* vs) {

    for (int i=0; i<max_num_masses; i++) {
        masses[i] = NULL;
    }
    for (int i=0; i<max_num_springs; i++) {
        springs[i] = NULL;
    }

    int masses_per_dim = ipow((VOX_SPACE_MAX_DEPTH+1)*2, 2);

    dfs_init_masses(vs,
                    0, 
                    masses, 
                    masses_per_dim, 
                    mass_count);

    init_springs(springs, masses, spring_count, max_num_masses, masses_per_dim);
    assert(*spring_count <= max_num_springs);

    for (int i=0; i<vs->num_voxels; i++) {


        for (int j=0; j<MAX_MASSES_PER_VOXEL; j++) {


        }
    }
}

int get_total_possible_masses(const int depth) {

    int dimension = (depth + 1) * 2;
    return ipow(dimension, 3);
}

int get_total_possible_springs(const int depth) {

    int dimension = (depth + 1) * 2;

    int inner_springs = 4 * total_voxels_at_depth(depth);
    int face_springs = 2 * dimension * ipow(dimension-1, 2);
    int edge_springs = 3 * (dimension - 1) * ipow(dimension, 2);

    return inner_springs + face_springs + edge_springs;
}

void dfs_init_masses(Voxel_space* vs,
                    const int idx, 
                    Mass** masses, 
                    const int masses_per_dim, 
                    int* mass_count) {

    if (1 == vs->tree[idx].exists) {

        for (int x=0; x<2; x++) {
            for (int y=0; y<2; y++) {
                for (int z=0; z<2; z++) {
                    
                    int mass_idx = vs->tree[idx].pos[0]+HEIGHT_OFFSET+x
                        + masses_per_dim*(vs->tree[idx].pos[1]+HEIGHT_OFFSET+y)
                        + ipow(masses_per_dim, 2)
                            *(vs->tree[idx].pos[2]+HEIGHT_OFFSET+z);

                    if (masses[mass_idx] == NULL) {

						(*mass_count)++;

            			masses[mass_idx] = malloc(sizeof(Mass));
            			masses[mass_idx]->m = MASS_M;
            			masses[mass_idx]->pos[0] = 
                            (vs->tree[idx].pos[0]+HEIGHT_OFFSET+x+0.0) * L0_SIDE;
						masses[mass_idx]->pos[0] = 
                            (vs->tree[idx].pos[1]+HEIGHT_OFFSET+y+0.0) * L0_SIDE;
						masses[mass_idx]->pos[0] = 
                            (vs->tree[idx].pos[2]+HEIGHT_OFFSET+z+0.0) * L0_SIDE;

            			for (int i=0; i<3; i++) {
                			masses[mass_idx]->vel[i] = 0.0f;
                			masses[mass_idx]->acc[i] = 0.0f;
            			}

                        masses[mass_idx]->material = vs->tree[idx].material;
        			}


                }  // end z
            }  // end y
        }  // end x

        voxel_type current_type = vs->tree[idx].type;

        if (ROOT == current_type) {
            for (int i=1; i<27; i++) {
                dfs_init_masses(vs, 
                                i, 
                                masses, 
                                masses_per_dim, 
                                mass_count);
            }
        }

        else if (MIDDLE == current_type) {
            dfs_init_masses(vs, 
                            get_child_index_of_m(idx), 
                            masses, 
                            masses_per_dim, 
                            mass_count);
        }

        else if (EDGE == current_type) {
            for (int i=0; i<2; i++) {
                dfs_init_masses(vs, 
                                get_child_index_of_e(idx, MIDDLE, i), 
                                masses, 
                                masses_per_dim, 
                                mass_count);
            }
              
            dfs_init_masses(vs, 
                            get_child_index_of_e(idx, EDGE, 0), 
                            masses, 
                            masses_per_dim, 
                            mass_count);
        }

        else if (CORNER == current_type) {
            for (int i=0; i<3; i++) {
                dfs_init_masses(vs, 
                                get_child_index_of_c(idx, MIDDLE, i), 
                                masses, 
                                masses_per_dim, 
                                mass_count);
            }

            for (int i=0; i<3; i++) {
                dfs_init_masses(vs, 
                                get_child_index_of_c(idx, EDGE, i), 
                                masses, 
                                masses_per_dim, 
                                mass_count);
            }

            dfs_init_masses(vs, 
                            get_child_index_of_c(idx, CORNER, 0), 
                            masses, 
                            masses_per_dim, 
                            mass_count);
        }

    }  // end if v->exists 
}

void init_springs(Spring** springs, 
                  Mass** masses, 
                  int* num_springs, 
                  const int max_num_masses,
                  const int masses_per_dim) {

    *num_springs = 0;

    for (int i=0; i<max_num_masses; i++) {
        if (NULL != masses[i]) {
            for (int x=0; x<2; x++) {
                for (int y=0; y<2; y++) {
                    for (int z=0; z<2; z++) {

                        // dont connect a mass to itself
                        if (x != 0 || y != 0 || z != 0) {
                            
                            int neighbor_idx = i 
                                + x 
                                + y * masses_per_dim
                                + z * ipow(masses_per_dim, 2);

                            if (neighbor_idx < max_num_masses
                                && NULL != masses[neighbor_idx]) {

                                springs[*num_springs] = malloc(sizeof(Spring));

                                springs[*num_springs]->m1 = i;
                                springs[*num_springs]->m2 = neighbor_idx;

                                float avg_k = 
                                    (material_to_k_map[masses[i]->material]
                        + material_to_k_map[masses[neighbor_idx]->material])/2.0f;

                                springs[*num_springs]->k = avg_k;
                                    
                                int dimension_diff = x + y + z;
                                float rest_length = length_map[dimension_diff];
                                springs[*num_springs]->l0 = rest_length;
                                springs[*num_springs]->a = rest_length; 

                                float avg_b;
                                float b1 = 
                                    get_b_from_mat(masses[i]->material,
                                                   springs[*num_springs]->l0);
                                float b2 = 
                                    get_b_from_mat(masses[neighbor_idx]->material,
                                                   springs[*num_springs]->l0);
                                avg_b = (b1 + b2)/2.0f;
            
                                springs[*num_springs]->b = avg_b;

                                float avg_c;
                                float c1 = 
                                    material_to_c_map[masses[i]->material];
                                float c2 = 
                                    material_to_c_map[
                                        masses[neighbor_idx]->material ];
                                avg_c = (c1 + c2)/2.0f;

                                springs[*num_springs]->c = avg_c;

                                (*num_springs)++;
                            }
                        }
                    }
                }
            }
        }
    }
}
