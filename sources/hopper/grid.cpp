#include <stdlib.h>
#include <string.h>
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
    // Idea is to see the array as the grid, and store particles into it
    // at their positions, converting to integers is nessecary (right?), and
    // should only slightly affect their real position so a few too many 
    // might be included.
    
    // I fear the whole double positioning might have eluded me though 
    // and this conversion wont work. Needs to be tested i guess.
    int coordinate = (int)particle->x * size + (int)particle->y;

    // Several particles may resolve to the same grid positions, so we save a 
    // grid position as a linked list of particles.
    linkedlist *gridPosParticle = (linkedlist *) malloc(sizeof(linkedlist));
    gridPosParticle->data = particle;

    // In case one already exist at the location we point to that (the order
    // they show up in is irellevant). IE we insert at beginning of the linked list.
    gridPosParticle->next = grid[coordinate];
    grid[coordinate] = gridPosParticle;
}

void grid_remove(particle_t *particle) {
    int coordinate = (int)particle->x * (int)size + (int)particle->y;

    linkedlist *listElem = grid[coordinate];
    linkedlist *prevElem;

    while(listElem && (listElem->data != particle)) {
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

