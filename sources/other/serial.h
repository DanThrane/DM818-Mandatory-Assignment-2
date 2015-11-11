#ifndef DM818_SERIAL_SERIAL_H
#define DM818_SERIAL_SERIAL_H

#include <stdio.h>
#include "common.h"

int runSerialImplementation(int argc, char **argv);
int runSerialWithParticles(FILE *fsave, int n, particle_t *particles);

#endif //DM818_SERIAL_SERIAL_H
