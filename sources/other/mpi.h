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
    particle_t *particles;
    int coordinateStart;
} GhostZone;

#endif //DM818_SERIAL_MPI_H
