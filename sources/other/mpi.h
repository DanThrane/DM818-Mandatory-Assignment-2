#ifndef DM818_SERIAL_MPI_H
#define DM818_SERIAL_MPI_H

#include <stdio.h>
#include "common.h"

int printHelp();

void initMPI(int &argc, char **&argv);

void initializeAtRoot(particle_t *particles, particle_t *particlesToSend, int *sendCount, int *sendDisplacement);

void validateScatter(particle_t *particlesToSend, particle_t *localParticles, int localCount);

void initSystem();

void validateRootInitialization(int currentProcessor, int counter, particle_t *particlesToSend);

void combineResult(FILE *fsave);

void wait_for_debugger();

typedef struct {
    int particleCount;
    particle_t *particles;
    int coordinateStart;
} GhostZone;

#ifdef DEBUG
#define VALIDATE_GHOST_ZONE(zone) {\
    assert(zone.particleCount >= 0);\
    for (int i = 0; i < zone.particleCount; i++) {\
        int coordinate = get_particle_coordinate(&zone.particles[i]);\
        assert(coordinate >= zone.coordinateStart && coordinate < zone.coordinateStart + gridColumns);\
    }\
}
#define WHEN_DEBUGGING(expr) expr
#define CASSERT(expr, message, ...) if (!(expr)) {\
    printf("[%d] Assertion failed at %s:%d. Message: " message "\n", rank, __FILE__, __LINE__, __VA_ARGS__);\
    assert(expr);\
}
#endif
#ifndef DEBUG
#define VALIDATE_GHOST_ZONE(zone)
#define WHEN_DEBUGGING(expr)
#define CASSERT(expr, message, ...)
#endif

#endif //DM818_SERIAL_MPI_H
