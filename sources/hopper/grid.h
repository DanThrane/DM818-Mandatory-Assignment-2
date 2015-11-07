#include "common.h"

struct linkedlist {
	linkedlist *next;
	particle_t *data;
} linkedlist;

void grid_init(double size);
void grid_add(particle_t *particle);
void grid_remove(particle_t *particle);
linkedlist* grid_get_collisions(particle_t *particle);
