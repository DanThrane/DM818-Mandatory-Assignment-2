#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include "grid.h"

linkedlist **grid;
double size;

void grid_init(double size) {
    size = size;
    grid = (linkedlist**) malloc(sizeof(linkedlist*) * (int)size**2);    
    memset(grid, 0, sizeof(linkedlist*) * (int)size**2);
}

void grid_add(particle_t *particle) {
    // Idea is to save the coordinate into the array as integers,
    // this should be correct enough to still find all the closest paritcles.
    // In the worst case this might just include a few too many?
    int coordinate = (int)particle->x * size + (int)particle->y;

    // Several particles may resolve to the same grid positions, so we save a 
    // grid position as a linked list.
    linkedlist *gridPosParticle = (linkedlist *) malloc(sizeof(linkedlist));
    gridPosParticle->value = particle;

    // In case one already exist at the location we point to that (the order
    // they show up in is irellevant).
    gridPosParticle->next = grid[coordinate];
    grid[coordinate] = gridPosParticle;
}

void grid_remove(particle_t *particle) {
    int coordinate = (int)particle->x * (int)size + (int)particle->y;

    linkedlist *listElem = grid[coordinate];
    linkedlist *prevElem;

    while(listElem && (listElem->value != particle)) {
        prevElem = listElem;
        listElem = listElem->next;
    }

    if (listElem) {
        prevElem->next = listElem->next;
        free(listElem);
    }
}

linkedlist* grid_get_collisions(particle_t *particle) {
    int coordinate = (int)particle->x * (int)size + (int)particle->y;
    linkedlist *listElem = grid[coordinate];
    return listElem;
}

