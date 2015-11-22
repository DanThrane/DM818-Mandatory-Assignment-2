#ifndef DM818_SERIAL_SERIAL_H
#define DM818_SERIAL_SERIAL_H

#include <stdio.h>
#include "common.h"

int runSerialImplementation(int argc, char **argv);
int runSerialWithParticles(FILE *fsave, int n, particle_t *particles);

#define DECLARE_TIMED_ZONE(name) double timer ## name = read_timer();
#define START_TIMED_ZONE(name) timer ## name = read_timer()
#define BEGIN_TIMED_ZONE(name) double timer ## name = read_timer()
#define END_TIMED_ZONE(name) double current ## name = 0; \
    auto result ## name = profiling.find(#name);\
    if (result ## name != profiling.end()) current ## name = result ## name->second;\
    double diff ## name = read_timer() - timer ## name; \
    profiling.insert({#name, current ## name + diff ## name})

#endif //DM818_SERIAL_SERIAL_H
