#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "grid.h"
#include "serial.h"

// benchmarking program
// note: i don't expect this to compile.
int runSerialImplementation(int argc, char **argv) {
    if (find_option(argc, argv, "-h") >= 0) {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);

    char *savename = read_string(argc, argv, "-o", NULL);

    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    particle_t *particles = (particle_t *) malloc(n * sizeof(particle_t));
    set_size(n);
    init_particles(n, particles);

    return runSerialWithParticles(fsave, n, particles);
}

int runSerialWithParticles(FILE *fsave, int n, particle_t *particles) {
    // new additions for linear runstime
    grid_init(sqrt(0.0005 * n) / 0.01);
    for (int i = 0; i < n; ++i) {
        grid_add(&particles[i]);
    }

    //  simulate a number of time steps
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {

        printf("NSTEP = %i\n", step);

        //  compute forces
        for (int i = 0; i < n; i++) {
            particles[i].ax = particles[i].ay = 0;

            // traverse included neighbors
            for (int offsetX = -1; offsetX <= 1; offsetX++) {
                for (int offsetY = -1; offsetY <= 1; offsetY++) {
                    const std::vector<particle_t *> &cell =
                            grid_get_collisions_at_loc((particles[i].x) + offsetX, (particles[i].y) + offsetY);

                    for (auto particle : cell) {
                        apply_force(particles[i], *particle);
                    }
                }
            }
        }

        //  move particles
        //  and update their position in the grid
        for (int i = 0; i < n; i++) {
            grid_remove(&particles[i]);
            move(particles[i]);
            grid_add(&particles[i]);
        }

        //  save if necessary
        if (fsave && (step % SAVEFREQ) == 0)
            save(fsave, n, particles);
    }
    simulation_time = read_timer() - simulation_time;

    printf("n = %d, simulation time = %g seconds\n", n, simulation_time);

    free(particles);
    if (fsave)
        fclose(fsave);

    return 0;
}
