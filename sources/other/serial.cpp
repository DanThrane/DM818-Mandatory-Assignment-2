#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <assert.h>

#include "common.h"
#include "grid.h"
#include "serial.h"

// benchmarking program
// note: i don't expect this to compile.
int runSerialImplementation(int argc, char **argv) {
#ifdef DEBUG
    printf("This is a debug build - This will have a large impact on the performance!\n");
#endif

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

std::vector<particle_t *> find_all_colliding_particles(int n, particle_t *particles, const particle_t &particle) {
    std::vector<particle_t *> result;
    for (int i = 0; i < n; i++) {
        particle_t *neighbor = &particles[i];

        double dx = neighbor->x - particle.x;
        double dy = neighbor->y - particle.y;
        double r2 = dx * dx + dy * dy;
        if (r2 <= cutoff * cutoff) {
            result.push_back(neighbor);
        }
    }
    return result;
}
int maxRows = 0;

int runSerialWithParticles(FILE *fsave, int n, particle_t *particles) {
    // new additions for linear runstime
    grid_init(sqrt(0.0005 * n) / 0.01);
    for (int i = 0; i < n; ++i) {
        grid_add(&particles[i]);
    }

    //  simulate a number of time steps
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {
#ifdef DEBUG
        validate_grid(n, __FILE__, __LINE__);
#endif
        printf("NSTEP = %i\n", step);

        //  compute forces
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

#ifdef DEBUG
            auto all = find_all_colliding_particles(n, particles, particles[i]);
            for (auto particle : all) {
                if (std::find(all_particles_visited.begin(), all_particles_visited.end(),
                              particle) == all_particles_visited.end()) {
                    printf("expected cell (%lf,%lf) to collide with (%lf,%lf) in iteration %i, it did not.\n",
                           particle->x, particle->y, particles[i].x, particles[i].y, step);
                    assert(false);
                }
            }
#endif
        }

        //  move particles
        //  and update their position in the grid
        for (int i = 0; i < n; i++) {
            int coordinate = get_particle_coordinate(&particles[i]);
            grid_remove(&particles[i]);
            move(particles[i]);
            grid_add(&particles[i]);
            int newCoordinate = get_particle_coordinate(&particles[i]);
            if (abs(coordinate - newCoordinate) > gridColumns) {
                int i1 = abs(coordinate - newCoordinate) / gridColumns;
//                printf("Moved %d cells (From %d to %d [%d rows])\n", abs(coordinate - newCoordinate), coordinate,
//                       newCoordinate, i1);
                if (i1 > maxRows) maxRows = i1;
            }
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
    printf("Max rows: %d\n", maxRows);

    return 0;
}
