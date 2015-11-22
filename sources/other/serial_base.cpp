#include <stdlib.h>
#include <tgmath.h>

#include "serial_base.h"
#include "grid.h"

int runSerialBaseImplementation(int argc, char **argv) {
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
    return runSerialBaseImplementationWithParticles(fsave, n, particles);
}

int maxRowsMoved = 0;

int runSerialBaseImplementationWithParticles(FILE *fsave, int n, particle_t *particles) {
    set_size(n);
    init_particles(n, particles);
    grid_init(sqrt(0.0005 * n) / 0.01);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {
        //
        //  compute forces
        //
        for (int i = 0; i < n; i++) {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++)
                apply_force(particles[i], particles[j]);
        }

        //
        //  move particles
        //
        for (int i = 0; i < n; i++) {
            int before = get_particle_coordinate(&particles[i]);
            move(particles[i]);
            int after = get_particle_coordinate(&particles[i]);
            int i1 = abs(before - after) / gridColumns;
            if (i1 > maxRowsMoved) maxRowsMoved = i1;
        }

        //
        //  save if necessary
        //
        if (fsave && (step % SAVEFREQ) == 0)
            save(fsave, n, particles);
    }
    simulation_time = read_timer() - simulation_time;

    printf("n = %d, simulation time = %g seconds\n", n, simulation_time);

    printf("Max rows: %d\n", maxRowsMoved);
    free(particles);
    if (fsave)
        fclose(fsave);

    return 0;
}
