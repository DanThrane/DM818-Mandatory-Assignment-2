#include <vector>

#include "common.h"

// Lets pretend this isn't C++

extern int gridColumns;

void gridInit(double size);

void gridAdd(particle_t *particle);

void gridRemove(particle_t *particle);

void gridValidate(int particle_count, const char *const file, int line);

std::vector<particle_t *> gridGetCollisionsAtNeighbor(particle_t *particle, int offsetX, int offsetY);

int gridGetParticleCoordinate(const particle_t *particle);

std::vector<particle_t *> gridGetAt(int coordinate);

void gridReset();

void gridPurge(int startInclusive, int endExclusive);

void gridValidateParticlesOnlyWithinZone(int startInclusive, int endInclusive);
