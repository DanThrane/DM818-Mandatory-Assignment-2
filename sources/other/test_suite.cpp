#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "common.h"
#include "serial_base.h"
#include "serial.h"

int main(int argc, char **argv) {
    if (find_option(argc, argv, "-h") >= 0) {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 100);

    FILE *fBase = fopen("base.txt", "w");
    FILE *fNew = fopen("new.txt", "w");

    particle_t *particles = (particle_t *) malloc(n * sizeof(particle_t));
    particle_t *particles2 = (particle_t *) malloc(n * sizeof(particle_t));

    set_size(n);
    init_particles(n, particles);
    memcpy(particles2, particles, n * sizeof(particle_t));

    runSerialBaseImplementationWithParticles(fBase, n, particles);
    reset();
//    runSerialBaseImplementationWithParticles(fNew, n, particles2);
//    runSerialWithParticles(fNew, n, particles2);
    return 0;
}
