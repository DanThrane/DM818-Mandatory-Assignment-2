#include "common.h"

struct linkedlist {
	linkedlist *next;
	particle_t *data;
} linkedlist;

void grid_init(int size);
void grid_add(particle_t *particle);
void grid_remove(particle_t *particle);
linkedlist_t* grid_get_collisions(particle_t *particle);
