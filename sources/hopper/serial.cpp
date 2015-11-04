#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "grid.h"

// benchmarking program
// note: i don't expect this to compile.
int main(int argc, char **argv) {
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

    // new additions for linear runstime
	double size = size = sqrt(0.0005 * n);
	grid_init(size);
	for (int i = 0; i < n; ++i) {
		grid_add(particles[i]);
	}

    //  simulate a number of time steps
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {
        
        //  compute forces
        for (int i = 0; i < n; i++) {
            particles[i].ax = particles[i].ay = 0;
            linkedlist *collisions = grid_get_collisions(particles[i]);
            while (collisions)
                // Currently this will ONLY calculate direct collisions,
                // need to be expanded to all grid locations around this too.
                apply_force(particles[i], collisions->data);
                collitions = collistions->next;
        }

        //  move particles
        //  todo: might have to remove/add particle to correct grid here.
        for (int i = 0; i < n; i++)
            move(particles[i]);

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
