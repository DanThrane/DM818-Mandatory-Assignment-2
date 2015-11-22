#ifndef DM818_SERIAL_SERIAL_BASE_H
#define DM818_SERIAL_SERIAL_BASE_H

#include <stdio.h>
#include "common.h"

int runSerialBaseImplementation(int argc, char **argv);
int runSerialBaseImplementationWithParticles(FILE *fsave, int n, particle_t *particles);

#endif //DM818_SERIAL_SERIAL_BASE_H
