/*
 * physics.c
 *
 */

#include <string.h>
#include <time.h>

#include "obj_parser.hpp"
#include "physics.h"

void init_masses_and_springs_from_voxel_space(Mass** masses, 
                                              const int max_num_masses,
                                              Spring** springs,
                                              const int max_num_springs,
                                              int* mass_count,
                                              int* spring_count,
                                              Voxel_space* vs,
                                              const float start_height) {

    for (int i=0; i<max_num_masses; i++) {
        masses[i] = NULL;
    }
    for (int i=0; i<max_num_springs; i++) {
        springs[i] = NULL;
    }

//    printf("max_num_masses: %d\n", max_num_masses);
//    printf("max_num_springs: %d\n", max_num_springs);

    int masses_per_dim = (VOX_SPACE_MAX_DEPTH+1)*2;

    dfs_init_masses(vs,
                    0, 
                    masses, 
                    masses_per_dim, 
                    mass_count,
                    start_height);


    init_springs(springs, masses, spring_count, max_num_masses, masses_per_dim);
//    printf("mass count: %d\n", *mass_count);
//    printf("spring count: %d\n", *spring_count);
    assert(*spring_count <= max_num_springs);

}

int get_total_possible_masses(const int depth) {

    int dimension = (depth + 1) * 2;
    return ipow(dimension, 3);
}

int get_total_possible_springs(const int depth) {

    int dimension = (depth + 1) * 2;

    int inner_crossing_springs = 4 * total_voxels_at_depth(depth);
    int horizontal_springs = dimension * (2*dimension*(dimension-1) 
                                          + 2*(dimension-1)*(dimension-1));
    int vertical_springs = (dimension-1) * ((dimension*dimension)
                                            + (2*dimension*(dimension-1)*2));
    /*
    int face_crossing_springs = 2 * dimension * (dimension-1)
            3 * dimension * ipow(dimension-1, 2);
    int edge_springs = 3 * (dimension - 1) * ipow(dimension, 2);*/

    return inner_crossing_springs + horizontal_springs + vertical_springs;
}

void dfs_init_masses(Voxel_space* vs,
                    const int idx, 
                    Mass** masses, 
                    const int masses_per_dim, 
                    int* mass_count,
                    const float start_height) {

    if (idx < total_voxels_at_depth(VOX_SPACE_MAX_DEPTH)
        && 1 == vs->tree[idx].exists) {

        for (int x=0; x<2; x++) {
            for (int y=0; y<2; y++) {
                for (int z=0; z<2; z++) {
                    
                    int mass_idx = vs->tree[idx].pos[0]+POS_OFFSET+x
                        + masses_per_dim*(vs->tree[idx].pos[1]+POS_OFFSET+y)
                        + ipow(masses_per_dim, 2)
                            *(vs->tree[idx].pos[2]+POS_OFFSET+z);

                    //printf("calculated mass idx: %d\n", mass_idx);

                    if (masses[mass_idx] == NULL) {

						(*mass_count)++;

            			masses[mass_idx] = (Mass*)malloc(sizeof(Mass));
                        CHECK_MALLOC_ERR(masses[mass_idx]);
            			masses[mass_idx]->m = MASS_M;
            			masses[mass_idx]->pos[0] = 
                            (vs->tree[idx].pos[0]+POS_OFFSET+x+0.0) * L0_SIDE;

						masses[mass_idx]->pos[1] = 
                            (vs->tree[idx].pos[1]+POS_OFFSET+y+0.0) * L0_SIDE;

						masses[mass_idx]->pos[2] = 
                            (vs->tree[idx].pos[2]+POS_OFFSET+z+0.0) * L0_SIDE
                            + start_height;

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
                                mass_count,
                                start_height);
            }
        }

        else if (MIDDLE == current_type) {
            dfs_init_masses(vs, 
                            get_child_index_of_m(idx), 
                            masses, 
                            masses_per_dim, 
                            mass_count,
                            start_height);
        }

        else if (EDGE == current_type) {
            for (int i=0; i<2; i++) {
                dfs_init_masses(vs, 
                                get_child_index_of_e(idx, MIDDLE, i), 
                                masses, 
                                masses_per_dim, 
                                mass_count,
                                start_height);
            }
              
            dfs_init_masses(vs, 
                            get_child_index_of_e(idx, EDGE, 0), 
                            masses, 
                            masses_per_dim, 
                            mass_count,
                            start_height);
        }

        else if (CORNER == current_type) {
            for (int i=0; i<3; i++) {
                dfs_init_masses(vs, 
                                get_child_index_of_c(idx, MIDDLE, i), 
                                masses, 
                                masses_per_dim, 
                                mass_count,
                                start_height);
            }

            for (int i=0; i<3; i++) {
                dfs_init_masses(vs, 
                                get_child_index_of_c(idx, EDGE, i), 
                                masses, 
                                masses_per_dim, 
                                mass_count,
                                start_height);
            }

            dfs_init_masses(vs, 
                            get_child_index_of_c(idx, CORNER, 0), 
                            masses, 
                            masses_per_dim, 
                            mass_count,
                            start_height);
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
                for (int y=-1; y<2; y++) {
                    for (int z=-1; z<2; z++) {

                        // dont connect a mass to itself
                        if ((x == 1) 
                            || (x == 0 && y == 1) 
                            || (x == 0 && y == 0 && z == 1)) {
                            
                            int neighbor_idx = i 
                                + x 
                                + y * masses_per_dim
                                + z * ipow(masses_per_dim, 2);

                            if (neighbor_idx >= 0
                                && neighbor_idx < max_num_masses
                                && NULL != masses[neighbor_idx]) {

                                float c2, c1, b2, b1, rest_length, avg_k;
                                int dimension_diff;

                                // this prevents us from accidentally wrapping 
                                // around:
                                for (int j=0; j<3; j++) {

                                    if (fabsf(masses[i]->pos[j] 
                                              - masses[neighbor_idx]->pos[j])
                                        > (L0_SIDE + 0.0001f)) {

                                        goto loop_end; 
                                    }  // end if fabsf
                                }  // end loop


                                springs[*num_springs] = (Spring*)malloc(sizeof(Spring));
                                CHECK_MALLOC_ERR(springs[*num_springs]);

                                springs[*num_springs]->m1 = i;
                                springs[*num_springs]->m2 = neighbor_idx;

                                avg_k = 
                                    (material_to_k_map[masses[i]->material]
                        + material_to_k_map[masses[neighbor_idx]->material])/2.0f;

                                springs[*num_springs]->k = avg_k;
                                    
                                dimension_diff = x + abs(y) + abs(z);
                                rest_length = length_map[dimension_diff];
                                assert(rest_length > 0.0f);
                                springs[*num_springs]->l0 = rest_length;
                                springs[*num_springs]->a = rest_length; 

                                float avg_b;
                                b1 = 
                                    get_b_from_mat(masses[i]->material,
                                                   springs[*num_springs]->l0);
                                b2 = 
                                    get_b_from_mat(masses[neighbor_idx]->material,
                                                   springs[*num_springs]->l0);
                                avg_b = (b1 + b2)/2.0f;

                                springs[*num_springs]->b = avg_b;

                                float avg_c;
                                c1 = 
                                    material_to_c_map[masses[i]->material];
                                c2 = 
                                    material_to_c_map[
                                        masses[neighbor_idx]->material ];
                                avg_c = (c1 + c2)/2.0f;

                                springs[*num_springs]->c = avg_c;

                                (*num_springs)++;

loop_end: ;
                            }  // end if neighbor_idx >=0
                        }  // end if x == 1
                    }  // end z loop
                }  // end y loop
            }  // end x loop
        }  // end if not null
    }  // end for each mass
}  // end function


/*
 * runs the physics simulation for each member of the population
 */
void simulate_population_cpu(Voxel_space** population, 
                             const int pop_size,
                             const float start_height) {

    // initialize masses and springs for each member of the population
    int max_masses_per_indiv = get_total_possible_masses(VOX_SPACE_MAX_DEPTH);
    int max_springs_per_indiv = get_total_possible_springs(VOX_SPACE_MAX_DEPTH);

    Mass*** pop_masses = (Mass***)malloc(sizeof(Mass**) * pop_size);
    CHECK_MALLOC_ERR(pop_masses);

    Spring*** pop_springs = (Spring***)malloc(sizeof(Spring**) * pop_size);
    CHECK_MALLOC_ERR(pop_springs);

    int* pop_mass_counts = (int*)malloc(sizeof(int) * pop_size);
    CHECK_MALLOC_ERR(pop_mass_counts);

    int* pop_spring_counts = (int*)malloc(sizeof(int) * pop_size);
    CHECK_MALLOC_ERR(pop_spring_counts);

    for (int i=0; i<pop_size; i++) {

        // init fitness to 0
        population[i]->fitness = 0.0f;

        Mass** indiv_masses = (Mass**)malloc(sizeof(Mass*) * max_masses_per_indiv);
        CHECK_MALLOC_ERR(indiv_masses);
        pop_masses[i] = indiv_masses;

        Spring** indiv_springs = (Spring**)malloc(sizeof(Spring*) * max_springs_per_indiv);
        CHECK_MALLOC_ERR(indiv_springs);
        pop_springs[i] = indiv_springs;

        pop_mass_counts[i] = 0;
        pop_spring_counts[i] = 0;
        init_masses_and_springs_from_voxel_space(pop_masses[i],
                                                 max_masses_per_indiv,
                                                 pop_springs[i],
                                                 max_springs_per_indiv,
                                                 &(pop_mass_counts[i]),
                                                 &(pop_spring_counts[i]),
                                                 population[i],
                                                 start_height);
    }
    
    float* force_vectors = (float*)malloc(sizeof(float) * pop_size*max_masses_per_indiv*3);
    CHECK_MALLOC_ERR(force_vectors);

    float* centers_of_mass_i = (float*)malloc(sizeof(float)*pop_size*3);
    CHECK_MALLOC_ERR(centers_of_mass_i);

    float* centers_of_mass_f = (float*)malloc(sizeof(float)*pop_size*3);
    CHECK_MALLOC_ERR(centers_of_mass_f);


    // get initial centers of mass:
    for (int i=0; i<pop_size; i++) {
        calculate_center_of_mass(pop_masses[i], 
                                 max_masses_per_indiv, 
                                 &(centers_of_mass_i[3*i]));
    }
    
    // this sucks but so do I
    int masses_per_indiv = max_masses_per_indiv;
   
    float t = 0.0f;
    for (int sim_i=0; sim_i<NUM_OF_ITERATIONS; sim_i++) {

        // reset forces
        for (int f=0; f<(pop_size*max_masses_per_indiv*3); f++) {
            force_vectors[f] = 0.0f;
        }
            
        // this assumes each cube has the same maximum dimensions !!!
        // spring loop:
        for (int i=0; i<(max_springs_per_indiv*pop_size); i++) {

            int indiv_idx = i / max_springs_per_indiv;
            int spring_idx = i % max_springs_per_indiv;

            // make sure spring exists (sorry not tabbing everything else again):
            if (NULL != pop_springs[indiv_idx][spring_idx]) {
                
            // increment each spring by its breathing function
            pop_springs[indiv_idx][spring_idx]->l0 = 
                pop_springs[indiv_idx][spring_idx]->a
                + pop_springs[indiv_idx][spring_idx]->b 
                * sinf(OMEGA*t + pop_springs[indiv_idx][spring_idx]->c);
            
#ifdef DEBUG
            printf("indiv: %d\n", indiv_idx);
            for (int j=0; j<max_masses_per_indiv; j++) {
                    
                if (pop_masses[indiv_idx][j] != NULL) {
                    printf("mass %d:\tx: %f, y: %f, z: %f\n", j,
                        pop_masses[indiv_idx][j]->pos[0],
                        pop_masses[indiv_idx][j]->pos[1],
                        pop_masses[indiv_idx][j]->pos[2]);

                }   
            }
            printf("\n");
#endif

                
                        
            // apply spring forces to masses
            int m1 = pop_springs[indiv_idx][spring_idx]->m1;
            int m2 = pop_springs[indiv_idx][spring_idx]->m2;
            float x1 = pop_masses[indiv_idx][m1]->pos[0];
            float y1 = pop_masses[indiv_idx][m1]->pos[1];
            float z1 = pop_masses[indiv_idx][m1]->pos[2];
            float x2 = pop_masses[indiv_idx][m2]->pos[0];
            float y2 = pop_masses[indiv_idx][m2]->pos[1];
            float z2 = pop_masses[indiv_idx][m2]->pos[2];
            
            float stretched_len = dist3d(x2, x1, y2, y1, z2, z1);

#ifdef DEBUG
            printf("indiv: %d, m1: %d, m2: %d\n", indiv_idx, m1, m2);
            printf("x1: %f, x2: %f, y1: %f, y2: %f, z1: %f, z2: %f\n", x1, x2, y1, y2, z1, z2);
            printf("stretched len: %f\n", stretched_len);
#endif

            //assert(stretched_len > 0.001);

            float force_normalized = 
                pop_springs[indiv_idx][spring_idx]->k
                * (stretched_len - pop_springs[indiv_idx][spring_idx]->l0);
            
            // update force vectors:
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m1 + 0]
                += force_normalized*(x2 - x1)/stretched_len;
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m1 + 1]
                += force_normalized*(y2 - y1)/stretched_len;
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m1 + 2]
                += force_normalized*(z2 - z1)/stretched_len;

            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m2 + 0]
                -= force_normalized*(x2 - x1)/stretched_len;
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m2 + 1]
                -= force_normalized*(y2 - y1)/stretched_len;
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m2 + 2]
                -= force_normalized*(z2 - z1)/stretched_len;

            }  // end if != NULL
        }

        // mass loop:
        // apply force due to gravity, check if each mass is on or below 
        // the ground and apply appropriate force, then update position
        for (int i=0; i<(max_masses_per_indiv*pop_size); i++) {

            int indiv_idx = i / max_masses_per_indiv;
            int mass_idx = i % max_masses_per_indiv;

            // make sure mass exists:
            if (NULL != pop_masses[indiv_idx][mass_idx]) {

#ifdef DEBUG
            if (0 == sim_i) {
                printf("starting positions for indiv %d mass %d\n", indiv_idx, mass_idx);
                printf("x: %f\n", pop_masses[indiv_idx][mass_idx]->pos[0]);
                printf("y: %f\n", pop_masses[indiv_idx][mass_idx]->pos[1]);
                printf("z: %f\n", pop_masses[indiv_idx][mass_idx]->pos[2]);
                printf("force vecs before calcs for indiv %d mass %d\n", indiv_idx, mass_idx);
                printf("force idx: %d\n", (masses_per_indiv*3)*indiv_idx+3*mass_idx+0);
                printf("x: %f\n", force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+0]);
                printf("y: %f\n", force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+1]);
                printf("z: %f\n", force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+2]);
            }
#endif

            // gravitational force:
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*mass_idx + 2]
                -= pop_masses[indiv_idx][mass_idx]->m * G;

            // ground forces (normal and friction):
            if (pop_masses[indiv_idx][mass_idx]->pos[2] <= 0) {

                // friction
                if (force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+2]<0) {

                    float force_horiz = 
    sqrtf(powf(force_vectors[(masses_per_indiv*3)*indiv_idx + 3*mass_idx + 0], 2.) 
    + powf(force_vectors[(masses_per_indiv*3)*indiv_idx + 3*mass_idx + 1], 2.));

                    float cos_theta = 
                    force_vectors[(masses_per_indiv*3)*indiv_idx + 3*mass_idx + 0]
                    / force_horiz;

                    float sin_theta = 
                    force_vectors[(masses_per_indiv*3)*indiv_idx + 3*mass_idx + 1]
                    / force_horiz;

                    float f_y = 
                    force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx + 2];

                    if (force_horiz < (-1) * f_y * U_S) {
                    force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+0]=0;
                    force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+1]=0;
                    }
                    else {
                        float f_kinetic_friction = U_K * f_y;
                        force_horiz += f_kinetic_friction;

                        force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+0] 
                            = force_horiz * cos_theta;
                        force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+1] 
                            = force_horiz * sin_theta;
                    }

                }  // end friction

                // normal force:
                force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx + 2]
                    = K_GROUND * fabsf(pop_masses[indiv_idx][mass_idx]->pos[2]);
            
            }  // end ground forces

#ifdef DEBUG
            if (0 == sim_i) {
                printf("first round of force vecs for indiv %d mass %d\n", indiv_idx, mass_idx);
                printf("x: %f\n", force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+0]);
                printf("y: %f\n", force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+1]);
                printf("z: %f\n", force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+2]);
            }
#endif


            // update positions:
            for (int j=0; j<3; j++) {
                // acceleration:
                pop_masses[indiv_idx][mass_idx]->acc[j] = 
                    force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx + j]
                    / pop_masses[indiv_idx][mass_idx]->m;

                // velocity:
                pop_masses[indiv_idx][mass_idx]->vel[j] +=
                    pop_masses[indiv_idx][mass_idx]->acc[j] * DT;

                pop_masses[indiv_idx][mass_idx]->vel[j] *= V_DAMP_CONST;

                // position:
                pop_masses[indiv_idx][mass_idx]->pos[j] +=
                    pop_masses[indiv_idx][mass_idx]->vel[j] * DT;                 
            }

            // apply height penalty to fitness:
            float threshold = ((2 * VOX_SPACE_MAX_DEPTH + 1)*L0_SIDE)+L0_SIDE;
            if (pop_masses[indiv_idx][mass_idx]->pos[2] > threshold) {
                population[indiv_idx]->fitness -= 
                    (pop_masses[indiv_idx][mass_idx]->pos[2] - threshold);
                    //* DT; 
            }

            }  // end mass exists check

        }
 
        t += DT;
    }  

    // get final centers of mass and calculate fitness:
    for (int i=0; i<pop_size; i++) {
        calculate_center_of_mass(pop_masses[i], 
                                 max_masses_per_indiv, 
                                 &(centers_of_mass_f[3*i]));

        population[i]->simulated_dist = centers_of_mass_f[3*i+0]
                                        - centers_of_mass_i[3*i+0];
        population[i]->fitness += population[i]->simulated_dist;
    }


    // clean up:

    free(force_vectors);
    free(centers_of_mass_i);
    free(centers_of_mass_f);

    for (int i=0; i<pop_size; i++) {
        for (int j=0; j<max_masses_per_indiv; j++) {
            if (NULL != pop_masses[i][j]) {
                free(pop_masses[i][j]);
            }
        }
        free(pop_masses[i]);

        for (int j=0; j<max_springs_per_indiv; j++) {
            if (NULL != pop_springs[i][j]) {
                free(pop_springs[i][j]);
            }
        }
        free(pop_springs[i]);
    }
    free(pop_masses);
    free(pop_springs);

    free(pop_mass_counts);
    free(pop_spring_counts);
}

void calculate_center_of_mass(Mass** indiv_masses, 
                              const int max_num_masses, 
                              float* com) {

    float total_mass = 0.0f;
    float pos_x = 0.0f;
    float pos_y = 0.0f;
    float pos_z = 0.0f;
    
    for (int i=0; i<max_num_masses; i++) {
        if (NULL != indiv_masses[i]) {
            total_mass += indiv_masses[i]->m;
            pos_x += indiv_masses[i]->m * indiv_masses[i]->pos[0];
            pos_y += indiv_masses[i]->m * indiv_masses[i]->pos[1];
            pos_y += indiv_masses[i]->m * indiv_masses[i]->pos[2];
        }
    }

    *(com+0) = pos_x / total_mass;
    *(com+1) = pos_y / total_mass;
    *(com+2) = pos_z / total_mass;
}

void export_to_gl(Voxel_space* vs, const float start_height) {

    int max_masses_per_indiv = get_total_possible_masses(VOX_SPACE_MAX_DEPTH);
    int max_springs_per_indiv = get_total_possible_springs(VOX_SPACE_MAX_DEPTH);

    Mass** indiv_masses = (Mass**)malloc(sizeof(Mass*) * max_masses_per_indiv);
    CHECK_MALLOC_ERR(indiv_masses);

    Spring** indiv_springs = (Spring**)malloc(sizeof(Spring*) * max_springs_per_indiv);
    CHECK_MALLOC_ERR(indiv_springs);

    int mass_count = 0;
    int spring_count = 0;

    init_masses_and_springs_from_voxel_space(indiv_masses,
                                             max_masses_per_indiv,
                                             indiv_springs,
                                             max_springs_per_indiv,
                                             &mass_count,
                                             &spring_count,
                                             vs,
                                             start_height);

    write_obj(vs, indiv_masses, max_masses_per_indiv, OBJ_FILE);

    char text[32];
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    strftime(text, sizeof(text), "%d-%H-%M", t);

    char walkfilename[64];
    strcpy(walkfilename, "vs-walk-");
    strcat(walkfilename, text);
    strcat(walkfilename, ".txt");

    char bouncefilename[64];
    strcpy(bouncefilename, "vs-bounce-");
    strcat(bouncefilename, text);
    strcat(bouncefilename, ".txt");

    //simulate_gl(vs, DEFAULT_START_HEIGHT, walkfilename);
    //simulate_gl(vs, 1.0, bouncefilename);

    free(indiv_masses);
    free(indiv_springs);
}

void write_obj(Voxel_space* vs,
               Mass** masses, 
               const int max_num_masses, 
               const char* filename) {

    FILE *f_obj = fopen(filename, "w");
    if (NULL == f_obj) {
        perror("could not open obj file for writing");
        exit(2);
    }

    fprintf(f_obj, "g voxel_space\n");

    for (int i=0; i<max_num_masses; i++) {
        if (NULL != masses[i]) {
            fprintf(f_obj, "v %f %f %f\n",
                     masses[i]->pos[0],
                     masses[i]->pos[1],
                     masses[i]->pos[2]);
        }
        else {
            // also print a dummy value for the unused masses so that we can
            // use the same indexing later
            fprintf(f_obj, "v 100.0 100.0 100.0\n");
        }
    }

    // vertex textures (i'm using them for colors lmao)
    fprintf(f_obj, "vt 0.5 0.0\n");
    fprintf(f_obj, "vt 1.0 1.0\n");
    fprintf(f_obj, "vt 1.0 0.5\n");
    fprintf(f_obj, "vt 0.1 0.5\n");

    // vertex normals
    fprintf(f_obj, "vn 0.0 0.0 1.0\n");
    fprintf(f_obj, "vn 0.0 1.0 0.0\n");
    fprintf(f_obj, "vn 0.0 0.0 -1.0\n");
    fprintf(f_obj, "vn 0.0 -1.0 0.0\n");
    fprintf(f_obj, "vn 1.0 0.0 0.0\n");
    fprintf(f_obj, "vn -1.0 0.0 0.0\n");

    int masses_per_dim = (VOX_SPACE_MAX_DEPTH+1)*2;

    int current_voxel_masses[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
    int total_voxels = total_voxels_at_depth(VOX_SPACE_MAX_DEPTH);
    for (int i=0; i<total_voxels; i++) {
        if (1 == vs->tree[i].exists) {
            
            int idx = 0;
            for (int x=0; x<2; x++) {
            for (int y=0; y<2; y++) {
            for (int z=0; z<2; z++) {

                int mass_idx = vs->tree[i].pos[0]+POS_OFFSET+x
                               + masses_per_dim*(vs->tree[i].pos[1]+POS_OFFSET+y)
                               + ipow(masses_per_dim, 2)
                                   *(vs->tree[i].pos[2]+POS_OFFSET+z);

                // plus 1 because obj files are 1-indexed
                current_voxel_masses[idx] = mass_idx+1;
                idx++;
            }
            }
            }

            int mat_color = vs->tree[i].material + 1;  // plus 1 bc obj 1-index

            // left side:
            fprintf(f_obj, "f %d/%d/6 ", 
                    current_voxel_masses[0],
                    mat_color);
            fprintf(f_obj, "%d/%d/6 ",
                    current_voxel_masses[1],
                    mat_color);
            fprintf(f_obj, "%d/%d/6\n",
                    current_voxel_masses[2],
                    mat_color);
            fprintf(f_obj, "f %d/%d/6 ", 
                    current_voxel_masses[2],
                    mat_color);
            fprintf(f_obj, "%d/%d/6 ",
                    current_voxel_masses[1],
                    mat_color);
            fprintf(f_obj, "%d/%d/6\n",
                    current_voxel_masses[3],
                    mat_color);

            // front side:
            fprintf(f_obj, "f %d/%d/4 ", 
                    current_voxel_masses[0],
                    mat_color);
            fprintf(f_obj, "%d/%d/4 ",
                    current_voxel_masses[5],
                    mat_color);
            fprintf(f_obj, "%d/%d/4\n",
                    current_voxel_masses[1],
                    mat_color);
            fprintf(f_obj, "f %d/%d/4 ", 
                    current_voxel_masses[0],
                    mat_color);
            fprintf(f_obj, "%d/%d/4 ",
                    current_voxel_masses[4],
                    mat_color);
            fprintf(f_obj, "%d/%d/4\n",
                    current_voxel_masses[5],
                    mat_color);

            // right side:
            fprintf(f_obj, "f %d/%d/5 ", 
                    current_voxel_masses[4],
                    mat_color);
            fprintf(f_obj, "%d/%d/5 ",
                    current_voxel_masses[7],
                    mat_color);
            fprintf(f_obj, "%d/%d/5\n",
                    current_voxel_masses[5],
                    mat_color);
            fprintf(f_obj, "f %d/%d/5 ", 
                    current_voxel_masses[4],
                    mat_color);
            fprintf(f_obj, "%d/%d/5 ",
                    current_voxel_masses[6],
                    mat_color);
            fprintf(f_obj, "%d/%d/5\n",
                    current_voxel_masses[7],
                    mat_color);

            // bottom side:
            fprintf(f_obj, "f %d/%d/3 ", 
                    current_voxel_masses[0],
                    mat_color);
            fprintf(f_obj, "%d/%d/3 ",
                    current_voxel_masses[2],
                    mat_color);
            fprintf(f_obj, "%d/%d/3\n",
                    current_voxel_masses[6],
                    mat_color);
            fprintf(f_obj, "f %d/%d/3 ", 
                    current_voxel_masses[0],
                    mat_color);
            fprintf(f_obj, "%d/%d/3 ",
                    current_voxel_masses[6],
                    mat_color);
            fprintf(f_obj, "%d/%d/3\n",
                    current_voxel_masses[4],
                    mat_color);

            // back side:
            fprintf(f_obj, "f %d/%d/2 ", 
                    current_voxel_masses[6],
                    mat_color);
            fprintf(f_obj, "%d/%d/2 ",
                    current_voxel_masses[3],
                    mat_color);
            fprintf(f_obj, "%d/%d/2\n",
                    current_voxel_masses[7],
                    mat_color);
            fprintf(f_obj, "f %d/%d/2 ", 
                    current_voxel_masses[6],
                    mat_color);
            fprintf(f_obj, "%d/%d/2 ",
                    current_voxel_masses[2],
                    mat_color);
            fprintf(f_obj, "%d/%d/2\n",
                    current_voxel_masses[3],
                    mat_color);

            // top side:
            fprintf(f_obj, "f %d/%d/1 ", 
                    current_voxel_masses[3],
                    mat_color);
            fprintf(f_obj, "%d/%d/1 ",
                    current_voxel_masses[5],
                    mat_color);
            fprintf(f_obj, "%d/%d/1\n",
                    current_voxel_masses[7],
                    mat_color);
            fprintf(f_obj, "f %d/%d/1 ", 
                    current_voxel_masses[3],
                    mat_color);
            fprintf(f_obj, "%d/%d/1 ",
                    current_voxel_masses[1],
                    mat_color);
            fprintf(f_obj, "%d/%d/1\n",
                    current_voxel_masses[5],
                    mat_color);
        }
    }

    // close file
    fprintf(f_obj, "\n");
    fclose(f_obj);
}

void simulate_gl(Voxel_space* vs, const float start_height, char* outfile) {

    

    FILE *f_out = fopen(outfile, "w");
    if (NULL == f_out) {
        perror("could not open obj file for writing");
        exit(2);
    }



    // TODO: not make this just a copy of the other simulate loop, but it works
    // for now
    int pop_size = 1;

    int max_masses_per_indiv = get_total_possible_masses(VOX_SPACE_MAX_DEPTH);
    int max_springs_per_indiv = get_total_possible_springs(VOX_SPACE_MAX_DEPTH);

    Mass*** pop_masses = (Mass***)malloc(sizeof(Mass**) * pop_size);
    CHECK_MALLOC_ERR(pop_masses);

    Spring*** pop_springs = (Spring***)malloc(sizeof(Spring**) * pop_size);
    CHECK_MALLOC_ERR(pop_springs);

    int* pop_mass_counts = (int*)malloc(sizeof(int) * pop_size);
    CHECK_MALLOC_ERR(pop_mass_counts);

    int* pop_spring_counts = (int*)malloc(sizeof(int) * pop_size);
    CHECK_MALLOC_ERR(pop_spring_counts);

    for (int i=0; i<pop_size; i++) {
        Mass** indiv_masses = (Mass**)malloc(sizeof(Mass*) * max_masses_per_indiv);
        CHECK_MALLOC_ERR(indiv_masses);
        pop_masses[i] = indiv_masses;

        Spring** indiv_springs = (Spring**)malloc(sizeof(Spring*) * max_springs_per_indiv);
        CHECK_MALLOC_ERR(indiv_springs);
        pop_springs[i] = indiv_springs;

        pop_mass_counts[i] = 0;
        pop_spring_counts[i] = 0;
        init_masses_and_springs_from_voxel_space(pop_masses[i],
                                                 max_masses_per_indiv,
                                                 pop_springs[i],
                                                 max_springs_per_indiv,
                                                 &(pop_mass_counts[i]),
                                                 &(pop_spring_counts[i]),
                                                 vs,
                                                 start_height);
    }

    float* force_vectors = (float*)malloc(sizeof(float) * pop_size*max_masses_per_indiv*3);
    CHECK_MALLOC_ERR(force_vectors);

    float* centers_of_mass_i = (float*)malloc(sizeof(float)*pop_size*3);
    CHECK_MALLOC_ERR(centers_of_mass_i);

    float* centers_of_mass_f = (float*)malloc(sizeof(float)*pop_size*3);
    CHECK_MALLOC_ERR(centers_of_mass_f);


    // get initial centers of mass:
    for (int i=0; i<pop_size; i++) {
        calculate_center_of_mass(pop_masses[i],
                                 max_masses_per_indiv,
                                 &(centers_of_mass_i[3*i]));
    }

    // this sucks but so do I
    int masses_per_indiv = max_masses_per_indiv;

    float t = 0.0f;
    for (int sim_i=0; sim_i<NUM_OF_ITERATIONS*2; sim_i++) {

        // reset forces
        for (int f=0; f<(pop_size*max_masses_per_indiv*3); f++) {
            force_vectors[f] = 0.0f;
        }

        // this assumes each cube has the same maximum dimensions !!!
        // spring loop:
        for (int i=0; i<(max_springs_per_indiv*pop_size); i++) {

            int indiv_idx = i / max_springs_per_indiv;
            int spring_idx = i % max_springs_per_indiv;

            // make sure spring exists (sorry not tabbing everything else again):
            if (NULL != pop_springs[indiv_idx][spring_idx]) {

            // increment each spring by its breathing function
            pop_springs[indiv_idx][spring_idx]->l0 =
                pop_springs[indiv_idx][spring_idx]->a
                + pop_springs[indiv_idx][spring_idx]->b
                * sinf(OMEGA*t + pop_springs[indiv_idx][spring_idx]->c);

            // apply spring forces to masses
            int m1 = pop_springs[indiv_idx][spring_idx]->m1;
            int m2 = pop_springs[indiv_idx][spring_idx]->m2;
            float x1 = pop_masses[indiv_idx][m1]->pos[0];
            float y1 = pop_masses[indiv_idx][m1]->pos[1];
            float z1 = pop_masses[indiv_idx][m1]->pos[2];
            float x2 = pop_masses[indiv_idx][m2]->pos[0];
            float y2 = pop_masses[indiv_idx][m2]->pos[1];
            float z2 = pop_masses[indiv_idx][m2]->pos[2];

            float stretched_len = dist3d(x2, x1, y2, y1, z2, z1);

            //assert(stretched_len > 0.001);

            float force_normalized =
                pop_springs[indiv_idx][spring_idx]->k
                * (stretched_len - pop_springs[indiv_idx][spring_idx]->l0);

            // update force vectors:
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m1 + 0]
                += force_normalized*(x2 - x1)/stretched_len;
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m1 + 1]
                += force_normalized*(y2 - y1)/stretched_len;
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m1 + 2]
                += force_normalized*(z2 - z1)/stretched_len;

            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m2 + 0]
                -= force_normalized*(x2 - x1)/stretched_len;
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m2 + 1]
                -= force_normalized*(y2 - y1)/stretched_len;
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*m2 + 2]
                -= force_normalized*(z2 - z1)/stretched_len;

            }  // end if != NULL
        }

        // mass loop:
        // apply force due to gravity, check if each mass is on or below
        // the ground and apply appropriate force, then update position
        for (int i=0; i<(max_masses_per_indiv*pop_size); i++) {

            int indiv_idx = i / max_masses_per_indiv;
            int mass_idx = i % max_masses_per_indiv;

            // make sure mass exists:
            if (NULL != pop_masses[indiv_idx][mass_idx]) {

            // gravitational force:
            force_vectors[(masses_per_indiv*3)*indiv_idx + 3*mass_idx + 2]
                -= pop_masses[indiv_idx][mass_idx]->m * G;

            // ground forces (normal and friction):
            if (pop_masses[indiv_idx][mass_idx]->pos[2] <= 0) {

                // friction
                if (force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+2]<0) {

                    float force_horiz =
    sqrtf(powf(force_vectors[(masses_per_indiv*3)*indiv_idx + 3*mass_idx + 0], 2.)
    + powf(force_vectors[(masses_per_indiv*3)*indiv_idx + 3*mass_idx + 1], 2.));

                    float cos_theta =
                    force_vectors[(masses_per_indiv*3)*indiv_idx + 3*mass_idx + 0]
                    / force_horiz;

                    float sin_theta =
                    force_vectors[(masses_per_indiv*3)*indiv_idx + 3*mass_idx + 1]
                    / force_horiz;

                    float f_y =
                    force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx + 2];

                    if (force_horiz < (-1) * f_y * U_S) {
                    force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+0]=0;
                    force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+1]=0;
                    }
                    else {
                        float f_kinetic_friction = U_K * f_y;
                        force_horiz += f_kinetic_friction;

                        force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+0]
                            = force_horiz * cos_theta;
                        force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx+1]
                            = force_horiz * sin_theta;
                    }

                }  // end friction

                // normal force:
                force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx + 2]
                    = K_GROUND * fabsf(pop_masses[indiv_idx][mass_idx]->pos[2]);

            }  // end ground forces

            // update positions:
            for (int j=0; j<3; j++) {
                // acceleration:
                pop_masses[indiv_idx][mass_idx]->acc[j] =
                    force_vectors[(masses_per_indiv*3)*indiv_idx+3*mass_idx + j]
                    / pop_masses[indiv_idx][mass_idx]->m;

                // velocity:
                pop_masses[indiv_idx][mass_idx]->vel[j] +=
                    pop_masses[indiv_idx][mass_idx]->acc[j] * DT;

                pop_masses[indiv_idx][mass_idx]->vel[j] *= V_DAMP_CONST;

                // position:
                pop_masses[indiv_idx][mass_idx]->pos[j] +=
                    pop_masses[indiv_idx][mass_idx]->vel[j] * DT;
            }

            }  // end mass exists check

        }

        if (sim_i%10==0) {
            write_obj(vs, 
                      pop_masses[0], 
                      max_masses_per_indiv,
                      "temp.txt");
            float* vp = NULL;
            float* vt = NULL;
            float* vn = NULL;
            int point_count = 0;
            load_obj_file("temp.txt", vp, vt, vn, point_count); 
            for (int p=0; p<point_count*3; p+=3) {
                fprintf(f_out, "%f,", vp[p]);
                fprintf(f_out, "%f,", vp[p+1]);
                fprintf(f_out, "%f", vp[p+2]);

                if (p != point_count*3 - 3) {
                    fprintf(f_out, ",");
                }

            }
            fprintf(f_out, "\n");
        }

        t += DT;
    }

    // clean up:

    fclose(f_out);

    free(force_vectors);
    free(centers_of_mass_i);
    free(centers_of_mass_f);

    for (int i=0; i<pop_size; i++) {
        for (int j=0; j<max_masses_per_indiv; j++) {
            if (NULL != pop_masses[i][j]) {
                free(pop_masses[i][j]);
            }
        }
        free(pop_masses[i]);

        for (int j=0; j<max_springs_per_indiv; j++) {
            if (NULL != pop_springs[i][j]) {
                free(pop_springs[i][j]);
            }
        }
        free(pop_springs[i]);
    }
    free(pop_masses);
    free(pop_springs);

    free(pop_mass_counts);
    free(pop_spring_counts);

    // close files
}

