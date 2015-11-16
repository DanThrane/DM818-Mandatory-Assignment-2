#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "grid.h"
#include "common.h"

struct linkedlist **grid;
int sizeD;

void grid_init(double size) {
    sizeD = size+1;
    grid = (struct linkedlist**) malloc(sizeof(struct linkedlist*) * (int)sizeD * (int)sizeD);    
    memset(grid, 0, sizeof(struct linkedlist*) * (int)sizeD * (int)sizeD);
}

void grid_add(particle_t *particle) {
    // Idea is to see the array as the grid, and store particles into it
    // at their positions, converting to integers is nessecary (right?), and
    // should only slightly affect their real position so a few too many 
    // might be included.
    
    // I fear the whole double positioning might have eluded me though 
    // and this conversion wont work. Needs to be tested i guess.
    int coordinate = (int)(particle->x/0.01) * sizeD + (int)(particle->y/0.01);

    // Several particles may resolve to the same grid positions, so we save a 
    // grid position as a linked list of particles.
    struct linkedlist * gridPosParticle = (struct linkedlist *) malloc(sizeof(linkedlist));
    gridPosParticle->data = particle;

    // In case one already exist at the location we point to that (the order
    // they show up in is irellevant). IE we insert at beginning of the linked list.
    gridPosParticle->next = grid[coordinate];
    grid[coordinate] = gridPosParticle;
}

void grid_remove(particle_t *particle) {
    int coordinate = (int)(particle->x/0.01) * (int)sizeD + (int)(particle->y/0.01);

    struct linkedlist *listElem = grid[coordinate];
    struct linkedlist *prevElem = grid[coordinate];

    while(listElem && (listElem->data != particle)) {
        prevElem = listElem;
        listElem = listElem->next;
    }

    if (listElem) {
        prevElem->next = listElem->next;
        //free(listElem);
    }
}

struct linkedlist* grid_get_collisions(particle_t *particle) {
    int coordinate = (int)(particle->x/0.01) * (int)sizeD + (int)(particle->y/0.01);
    return grid[coordinate];
}

struct linkedlist* grid_get_collisions_at_loc(double x, double y) {
    // Handle incorrect positions:
    if (x/0.01 < 0 || y/0.01 < 0)
        return NULL;
    if (x/0.01 > sizeD-1 || y/0.01 > sizeD-1)
        return NULL;

    int coordinate = (int)(x/0.01) * (int)sizeD + (int)(y/0.01);
    return grid[coordinate];
}