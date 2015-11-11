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
	grid_init(sqrt(0.0005 * n)/0.01);
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
            linkedlist *collisions = grid_get_collisions(&particles[i]);
            while (collisions != 0) {
                // TODO: Currently this will ONLY calculate direct collisions,
                // need to be expanded to all grid locations / neighbors around 
		// TODO: Segfaults because while check dont seem to work?
		apply_force(particles[i], *(collisions->data));
                collisions = collisions->next;
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
