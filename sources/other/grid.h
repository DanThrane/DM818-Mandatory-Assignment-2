#include <vector>


#include "common.h"

extern int gridColumns;

void grid_init(double size);

void grid_add(particle_t *particle);

void grid_remove(particle_t *particle);

void validate_grid(int particle_count);

std::vector<particle_t *> grid_get_collisions(particle_t *particle);

std::vector<particle_t *> grid_get_collisions_at_neighbor(particle_t *particle, int offsetX, int offsetY);

int get_particle_coordinate(const particle_t *particle);

std::vector<particle_t *> grid_get_at(int coordinate);

void grid_reset();

void grid_purge(int startInclusive, int endExclusive);

void grid_disable_track();

void grid_enable_track();

void grid_track_insertions(int startInclusive, int endExclusive, std::vector<particle_t *> *consumer);