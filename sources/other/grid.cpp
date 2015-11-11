#include <vector>
#include <algorithm>
#include <assert.h>

#include "grid.h"

std::vector<std::vector<particle_t *> > grid;
int sizeD;

int get_particle_coordinate(const particle_t *particle);

void grid_init(double size) {
    sizeD = (int) (size + 1);
    grid.resize((unsigned long) (sizeD * sizeD));
}

void validate_grid(int particle_count) {
    int particles_in_system = 0;
    for (int i = 0; i < grid.size(); i++) {
        for (auto particle : grid[i]) {
            assert(get_particle_coordinate(particle) == i);
        }
        particles_in_system += grid[i].size();
    }
    assert(particles_in_system == particle_count);
}

void grid_add(particle_t *particle) {
    // Idea is to see the array as the grid, and store particles into it
    // at their positions, converting to integers is nessecary (right?), and
    // should only slightly affect their real position so a few too many 
    // might be included.

    // I fear the whole double positioning might have eluded me though 
    // and this conversion wont work. Needs to be tested i guess.
    int coordinate = get_particle_coordinate(particle);

    grid[coordinate].push_back(particle);
}

int get_particle_coordinate(const particle_t *particle) {
    int coordinate = (int) (particle->x / 0.01) * sizeD + (int) (particle->y / 0.01);
    return coordinate;
}

void grid_remove(particle_t *particle) {
    int coordinate = get_particle_coordinate(particle);

    std::vector<particle_t *> &cell = grid[coordinate];
    cell.erase(std::remove(cell.begin(), cell.end(), particle), cell.end());
}

std::vector<particle_t *> grid_get_collisions(particle_t *particle) {
    int coordinate = get_particle_coordinate(particle);
    return grid[coordinate];
}

std::vector<particle_t *> grid_get_collisions_at_loc(double x, double y) {
    // Handle incorrect positions:
    if (x / 0.01 < 0 || y / 0.01 < 0)
        return std::vector<particle_t *>();
    if (x / 0.01 > sizeD - 1 || y / 0.01 > sizeD - 1)
        return std::vector<particle_t *>();

    int coordinate = (int) (x / 0.01) * sizeD + (int) (y / 0.01);
    return grid[coordinate];
}
