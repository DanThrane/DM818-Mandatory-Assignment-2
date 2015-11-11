#include <vector>

#include "common.h"

void grid_init(double size);

void grid_add(particle_t *particle);

void grid_remove(particle_t *particle);

void validate_grid(int particle_count);

std::vector<particle_t *> grid_get_collisions(particle_t *particle);

std::vector<particle_t *> grid_get_collisions_at_loc(double x, double y);
