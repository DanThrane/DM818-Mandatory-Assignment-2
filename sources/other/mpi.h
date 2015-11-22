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

void waitForDebugger();

typedef struct {
    int particleCount;
    particle_t *particles;
    int coordinateStart;
} GhostZone;

typedef struct {
    int id;
    int step;
    double x;
    double y;
} SavedParticle;

#ifdef DEBUG
#define VALIDATE_GHOST_ZONE(zone) {\
    assert(zone.particleCount >= 0);\
    for (int i = 0; i < zone.particleCount; i++) {\
        int coordinate = get_particle_coordinate(&zone.particles[i]);\
        assert(coordinate >= zone.coordinateStart && coordinate < zone.coordinateStart + gridColumns);\
    }\
}
#define WHEN_DEBUGGING(expr) expr
#endif
#ifndef DEBUG
#define VALIDATE_GHOST_ZONE(zone)
#define WHEN_DEBUGGING(expr)
#endif

#endif //DM818_SERIAL_MPI_H
