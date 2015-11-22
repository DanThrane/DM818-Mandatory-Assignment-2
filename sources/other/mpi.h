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

#define DECLARE_TIMED_ZONE(name) double timer ## name = read_timer();
#define START_TIMED_ZONE(name) timer ## name = read_timer()
#define BEGIN_TIMED_ZONE(name) double timer ## name = read_timer()
#define END_TIMED_ZONE(name) double current ## name = 0; \
    auto result ## name = profiling.find(#name);\
    if (result ## name != profiling.end()) current ## name = result ## name->second;\
    double diff ## name = read_timer() - timer ## name; \
    profiling.insert({#name, current ## name + diff ## name})

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
