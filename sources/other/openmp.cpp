#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <algorithm>
#include "common.h"
#include "grid.h"

//
//  benchmarking program
//
int main(int argc, char **argv) {
    if (find_option(argc, argv, "-h") >= 0) {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);
    char *savename = read_string(argc, argv, "-o", NULL);

    FILE *fsave = savename ? fopen(savename, "w") : NULL;

    particle_t *particles = (particle_t *) malloc(n * sizeof(particle_t));
    set_size(n);
    init_particles(n, particles);

    // new additions for linear runtime
    // could pragma this, but its so simple there should be no need.
    grid_init(sqrt(0.0005 * n) / 0.01);
    for (int i = 0; i < n; ++i) {
        grid_add(&particles[i]);
    }

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();

    #pragma omp parallel
    for (int step = 0; step < 1000; step++) {
        printf("NSTEP = %i\n", step);

        #ifdef DEBUG
                validate_grid(n);
        #endif

        //
        //  compute all forces
        //
        #pragma omp for
        for (int i = 0; i < n; i++) {
            particles[i].ax = particles[i].ay = 0;

            // traverse included neighbors
            #ifdef DEBUG
                        std::vector<particle_t *> all_particles_visited;
            #endif
            for (int offsetX = -1; offsetX <= 1; offsetX++) {
                for (int offsetY = -1; offsetY <= 1; offsetY++) {
                    const std::vector<particle_t *> &cell =
                            grid_get_collisions_at_neighbor(&particles[i], offsetX, offsetY);

                    for (auto particle : cell) {
                        apply_force(particles[i], *particle);
                        #ifdef DEBUG
                                                all_particles_visited.push_back(particle);
                        #endif
                    }
                }
            }
        }

        //
        //  move particles
        //
        #pragma omp for
        for (int i = 0; i < n; i++) {
            grid_remove(&particles[i]);
            move(particles[i]);
            grid_add(&particles[i]);
        }

        //
        //  save if necessary
        //
        #pragma omp master
        if (fsave && (step % SAVEFREQ) == 0) {
            save(fsave, n, particles);
        }
    }
    simulation_time = read_timer() - simulation_time;

    printf("n = %d,\tsimulation time = %g seconds\n", n, simulation_time);

    free(particles);
    if (fsave)
        fclose(fsave);

    return 0;
}
