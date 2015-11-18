#include <vector>
#include <algorithm>
#include <assert.h>

#include "grid.h"

std::vector<std::vector<particle_t *> > grid;
int gridColumns;

int get_particle_coordinate(const particle_t *particle);

void grid_init(double size) {
    gridColumns = (int) (size + 1);
    grid.resize((unsigned long) (gridColumns * gridColumns));
}

void grid_reset() {
    grid.clear();
    grid.resize((unsigned long) (gridColumns * gridColumns));
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

int get_particle_coordinate(double x, double y) {
    int coordinate = (int) (x / 0.01) * gridColumns + (int) (y / 0.01);
    return coordinate;
}

int get_particle_coordinate(const particle_t *particle) {
    int coordinate = (int) (particle->x / 0.01) * gridColumns + (int) (particle->y / 0.01);
    return coordinate;
}

void grid_remove(particle_t *particle) {
    int coordinate = get_particle_coordinate(particle);

    std::vector<particle_t *> &cell = grid[coordinate];
    cell.erase(std::remove(cell.begin(), cell.end(), particle), cell.end());
}

void grid_purge(int startInclusive, int endExclusive) {
    for (int i = startInclusive; i < endExclusive; i++) {
        // TODO Release particles too?
        grid[i].clear();
    }
}

std::vector<particle_t *> grid_get_collisions(particle_t *particle) {
    int coordinate = get_particle_coordinate(particle);
    return grid[coordinate];
}

std::vector<particle_t *> grid_get_at(int coordinate) {
    return grid[coordinate];
}

std::vector<particle_t *> grid_get_collisions_at_neighbor(particle_t *particle, int offsetX, int offsetY) {
    int coordinate = get_particle_coordinate(particle);
    int x = (coordinate % gridColumns) + offsetX;
    int y = (coordinate / gridColumns) + offsetY;

    if (x < 0 || y < 0 || x >= gridColumns || y >= gridColumns) {
        return std::vector<particle_t *>();
    }

    return grid[coordinate + (offsetY * 23) + offsetX];
}
