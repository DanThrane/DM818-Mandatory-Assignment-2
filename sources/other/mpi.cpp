#include <mpi.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sstream>

#include "grid.h"
#include "mpi.h"

particle_t shitParticle1;
particle_t shitParticle2;

//
// Global state
//

/**
 * Output stream (for particle locations)
 */
std::ostringstream localOutput;

/**
 * The number of cells for each process
 */
int cellsPerProcess;

/**
 * The rank of this processor
 */
int rank;

/**
 * The maximum rank (this is the total number of processors in the system)
 */
int maxRank;

/**
 * The particle data type used by MPI
 */
MPI_Datatype particleType;

/**
 * The amount of particles that exist in the entire system
 */
int globalParticleCount;

//
// Local state
//

/**
 * The amount of particles owned by this processor
 */
int maxPosition;

/**
 * The particles owned by this processor
 */
particle_t *ownedParticles;

GhostZone ownedUpper;
GhostZone ownedLower;

GhostZone borrowedUpper;
GhostZone borrowedLower;

std::vector<particle_t> insertionsIntoUpperBorrowed;
std::vector<particle_t> insertionsIntoLowerBorrowed;

void validateNumberOfParticles() {
    int count;
    MPI_Reduce(&maxPosition, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        CASSERT(count == globalParticleCount, "Expected %d particles, but only found %d", globalParticleCount, count);
    }
}

void prepareGhostZoneForExchange(GhostZone &zone) {
    // Move the particles from the real grid into the buffer
    zone.particleCount = 0;
    for (int i = zone.coordinateStart; i < zone.coordinateStart + gridColumns; i++) {
        for (auto particle : grid_get_at(i)) {
            memcpy(&zone.particles[zone.particleCount], particle, sizeof(particle_t));
            zone.particleCount++;
#ifdef DEBUG
            int coordinate = get_particle_coordinate(particle);
            if (!(coordinate >= zone.coordinateStart && coordinate < zone.coordinateStart + gridColumns)) {
                printf("[%d] Found particle at %d, but not supposed to exceed zone starting at %d. "
                               "(gridCols = %d, diff = %d)\n", rank, coordinate, zone.coordinateStart, gridColumns,
                       coordinate - zone.coordinateStart);
                assert(false);
            }
#endif
        }
    }
}

void exchangeParticles(GhostZone &owned, GhostZone &borrowed, int multiplier) {
    // Communicate number of particles in update
    MPI_Sendrecv(&owned.particleCount, 1, MPI_INT, rank + (1 * multiplier), 0, &borrowed.particleCount, 1,
                 MPI_INT, rank + (1 * multiplier), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Exchange particles
    MPI_Sendrecv(owned.particles, owned.particleCount, particleType, rank + (1 * multiplier), 0,
                 borrowed.particles, borrowed.particleCount, particleType, rank + (1 * multiplier), 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    VALIDATE_GHOST_ZONE(owned);
    VALIDATE_GHOST_ZONE(borrowed);
}

particle_t *prepareInsertions(std::vector<particle_t> insertions) {
    int i = 0;
    particle_t *result = (particle_t *) malloc(sizeof(particle_t) * insertions.size());
    for (auto particle : insertions) {
        memcpy(&result[i], &particle, sizeof(particle_t));
        i++;
    }
    return result;
}

particle_t *exchangeInsertions(std::vector<particle_t> &insertions, int multiplier, int *outputCount) {
    int count;
    int sendingCount = (int) insertions.size();

    particle_t *prepared = prepareInsertions(insertions);
    particle_t *receiveBuffer;

    MPI_Sendrecv(&sendingCount, 1, MPI_INT, rank + (1 * multiplier), 0, &count, 1, MPI_INT, rank + (1 * multiplier), 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    receiveBuffer = (particle_t *) malloc(sizeof(particle_t) * count);
//    printf("[%d] Receiving %d inserted particles!\n", rank, count);
    MPI_Sendrecv(prepared, sendingCount, particleType, rank + (1 * multiplier), 0, receiveBuffer, count, particleType,
                 rank + (1 * multiplier), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    *outputCount = count;
    free(prepared);
    return receiveBuffer;
}

void updateGrid(GhostZone &zone, std::vector<particle_t> &localInsertions) {
    VALIDATE_GHOST_ZONE(zone);
    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));

    for (int i = 0; i < zone.particleCount; i++) {
        grid_add(&zone.particles[i]);
    }

    for (auto particle : localInsertions) {
        mempcpy(&zone.particles[zone.particleCount], &particle, sizeof(particle_t));
        grid_add(&zone.particles[zone.particleCount]);
        zone.particleCount++;
    }

    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
    VALIDATE_GHOST_ZONE(zone);
}

void mergeInsertedInLocallyOwned(int insertedCount, particle_t *inserted) {
    int start = maxPosition;

    // The incoming buffer is a temporary one. Move the particles to the owned buffer
    memcpy(&ownedParticles[maxPosition], inserted, sizeof(particle_t) * insertedCount);
    maxPosition += insertedCount;
    if (rank == 0) printf(" += %d\n", insertedCount);

    for (int i = start; i < start + insertedCount; i++) {
        grid_add(&ownedParticles[i]);
    }
}

void exchangeInformationAbove(particle_t **insertedLower, int *insertedLowerCount) {
    if (rank < maxRank - 1) {
        WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
        exchangeParticles(ownedUpper, borrowedUpper, 1);
        WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
        *insertedLower = exchangeInsertions(insertionsIntoUpperBorrowed, 1, insertedLowerCount);
        WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
    }
}

void exchangeInformationBelow(particle_t **insertedUpper, int *insertedUpperCount) {
    if (rank > 0) {
        WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
        exchangeParticles(ownedLower, borrowedLower, -1);
        WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
        *insertedUpper = exchangeInsertions(insertionsIntoLowerBorrowed, -1, insertedUpperCount);
        WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
    }
}

void exchangeInformationWithNeighborHood() {
    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
#ifdef DEBUG
    // Validate our local insertions, which are to be exchanged with neighbor
    int lowerStart = borrowedLower.coordinateStart;
    int lowerEnd = lowerStart + gridColumns;
    for (auto particle : insertionsIntoLowerBorrowed) {
        int coordinate = get_particle_coordinate(&particle);
        assert(coordinate >= lowerStart && coordinate < lowerEnd);
    }

    int upperStart = borrowedUpper.coordinateStart;
    int upperEnd = upperStart + gridColumns;
    for (auto particle : insertionsIntoUpperBorrowed) {
        int coordinate = get_particle_coordinate(&particle);
        assert(coordinate >= upperStart && coordinate < upperEnd);
    }
#endif
    int insertedIntoUpperOwnedCount = 0;
    particle_t *insertedIntoUpperOwned = NULL;
    int insertedIntoLowerOwnedCount = 0;
    particle_t *insertedIntoLowerOwned = NULL;

    // Prepare ghost zones such that we're ready to exchange them with our neighbors
    prepareGhostZoneForExchange(ownedUpper);
    prepareGhostZoneForExchange(ownedLower);

    // Purge the borrowed ghost zones, as we will receive new particles soon.
    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
    grid_purge(borrowedLower.coordinateStart, borrowedLower.coordinateStart + gridColumns);
    grid_purge(borrowedUpper.coordinateStart, borrowedUpper.coordinateStart + gridColumns);
    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));

    if (rank % 2 == 0) {
        // Even ranks communicate with processors above us first
        exchangeInformationAbove(&insertedIntoLowerOwned, &insertedIntoLowerOwnedCount);
        exchangeInformationBelow(&insertedIntoUpperOwned, &insertedIntoUpperOwnedCount);
    } else {
        // Even ranks communicate with processors below us first
        exchangeInformationBelow(&insertedIntoUpperOwned, &insertedIntoUpperOwnedCount);
        exchangeInformationAbove(&insertedIntoLowerOwned, &insertedIntoLowerOwnedCount);
    }

    // First insert the complete state changes from our neighbor merged with out own local insertions in these zones
    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
    updateGrid(borrowedUpper, insertionsIntoUpperBorrowed);
    updateGrid(borrowedLower, insertionsIntoLowerBorrowed);
    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));

    // Then we get the insertions that occurred into our owned zones from our neighbors
    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
    mergeInsertedInLocallyOwned(insertedIntoLowerOwnedCount, insertedIntoLowerOwned);
    mergeInsertedInLocallyOwned(insertedIntoUpperOwnedCount, insertedIntoUpperOwned);
    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));

    // Clear insertions from last iteration
    insertionsIntoLowerBorrowed.clear();
    insertionsIntoUpperBorrowed.clear();

    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));

    // Release inserted
    if (insertedIntoLowerOwned == NULL) free(insertedIntoLowerOwned);
    if (insertedIntoUpperOwned == NULL) free(insertedIntoUpperOwned);
}

int main(int argc, char **argv) {
    //  process command line parameters
    globalParticleCount = read_int(argc, argv, "-n", 1000);
    char *savename = read_string(argc, argv, "-o", NULL);
    FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : NULL;
    if (find_option(argc, argv, "-h") >= 0) { return printHelp(); }

    initMPI(argc, argv);
    initSystem();
    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));

//    if (rank == 0) wait_for_debugger();

    //  simulate a number of time steps
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {
        if (rank == 0) printf("Step = %d. MaxPosition = %d\n", step, maxPosition);
        MPI_Barrier(MPI_COMM_WORLD);
        WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
        exchangeInformationWithNeighborHood();
        WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
        WHEN_DEBUGGING(validateNumberOfParticles());

        for (int i = 0; i < gridColumns * gridColumns; i++) {
            if (i < borrowedLower.coordinateStart || i > borrowedUpper.coordinateStart + gridColumns) {
                if (!grid_get_at(i).empty()) {
                    printf("Not empty at %d (It has %d particles)\n", i, (int) grid_get_at(i).size());
                }
            }
        }

        //  compute forces
        for (int i = 0; i < maxPosition; i++) {
            ownedParticles[i].ax = ownedParticles[i].ay = 0;

            // traverse included neighbors
            for (int offsetX = -1; offsetX <= 1; offsetX++) {
                for (int offsetY = -1; offsetY <= 1; offsetY++) {
                    const std::vector<particle_t *> &cell =
                            grid_get_collisions_at_neighbor(&ownedParticles[i], offsetX, offsetY);

                    for (auto particle : cell) {
                        if (&ownedParticles[i] == particle) continue; //This will do very bad things otherwise!!!
                        memcpy(&shitParticle1, particle, sizeof(particle_t));
                        memcpy(&shitParticle2, &ownedParticles[i], sizeof(particle_t));

                        apply_force(ownedParticles[i], *particle);
                    }
                }
            }
        }

        // Move and update particles
        for (int i = 0; i < maxPosition; i++) {
#ifdef DEBUG
            int min = borrowedLower.coordinateStart;
            int max = borrowedUpper.coordinateStart + gridColumns;

            if (min < 0) min = cellsPerProcess * rank;
            if (max >= gridColumns * gridColumns) max = (gridColumns * gridColumns);

            validate_particles_within_sub_grid(min, max);

            particle_t debugParticleBefore;
#endif
            particle_t *particle = &ownedParticles[i];
#ifdef DEBUG
            memcpy(&debugParticleBefore, particle, sizeof(particle_t));
#endif
            int beforeStart = get_particle_coordinate(particle);
            grid_remove(particle);
            move(*particle);
            WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
            int coordinate = get_particle_coordinate(particle);
            if (coordinate < borrowedLower.coordinateStart + gridColumns ||
                coordinate >= borrowedUpper.coordinateStart) { // Out of bounds
                // Copy the particle and insert into a vector of insertions into a certain ghost zone
                particle_t copy;
                memcpy(&copy, particle, sizeof(particle_t));

                if (coordinate >= borrowedLower.coordinateStart &&
                    coordinate < borrowedLower.coordinateStart + gridColumns) {
                    insertionsIntoLowerBorrowed.push_back(copy);
                } else if (coordinate >= borrowedUpper.coordinateStart &&
                           coordinate < borrowedUpper.coordinateStart + gridColumns) {
                    insertionsIntoUpperBorrowed.push_back(copy);
                } else {
#ifdef DEBUG
                    WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
                    VALIDATE_GHOST_ZONE(borrowedLower);
                    VALIDATE_GHOST_ZONE(borrowedUpper);

                    printf("Coordinate: %d, before: %d, lower: %d upper: %d\n", coordinate, beforeStart,
                           borrowedLower.coordinateStart, borrowedUpper.coordinateStart);
                    printf("Particle before: p.x=%f, p.y=%f, p.ax=%f, p.ay=%f, p.vx=%f, p.vy=%f\n",
                           debugParticleBefore.x, debugParticleBefore.y, debugParticleBefore.ax, debugParticleBefore.ay,
                           debugParticleBefore.vx, debugParticleBefore.vy);
                    printf("Particle after: p.x=%f, p.y=%f, p.ax=%f, p.ay=%f, p.vx=%f, p.vy=%f\n",
                           particle->x, particle->y, particle->ax, particle->ay,
                           particle->vx, particle->vy);
                    assert(false);
#endif
                }

                // Take last element and move it in its place
                maxPosition--;
                printf("--\n");
                if (i != maxPosition) {
                    grid_remove(&ownedParticles[maxPosition]); // Make sure that the grid gets the updated pointer
                    memcpy(&ownedParticles[i], &ownedParticles[maxPosition], sizeof(particle_t));
                    grid_add(&ownedParticles[i]);
                }
            } else {
                grid_add(particle);
            }
            WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
        }

        //  save if necessary
        if (savename && step % SAVEFREQ == 0) {
            for (int i = 0; i < maxPosition; i++) {
                localOutput << ownedParticles[i].x << " " << ownedParticles[i].y << "\n";
            }
        }
        // Communicate to rank 0 with result
    }
    simulation_time = read_timer() - simulation_time;

    if (savename) {
        combineResult(fsave);
    }

    if (rank == 0) { printf("n = %d, simulation time = %g seconds\n", globalParticleCount, simulation_time); }

    if (rank == 0 && fsave) {
        fprintf(fsave, localOutput.str().c_str());
    }

    // TODO Release resources
    if (fsave) {
        fclose(fsave);
    }

    MPI_Finalize();
    return 0;
}

void combineResult(FILE *fsave) {
    int counts[maxRank];
    char *combinedOutput;
    int *displacements;

    int count = 0;
    std::string output = localOutput.str();
    int outputSize = (int) output.size();

    MPI_Gather(&outputSize, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        displacements = (int *) malloc(sizeof(int) * maxRank);
        for (int i = 0; i < maxRank; i++) {
            displacements[i] = count;
            count += counts[i];
        }
        combinedOutput = (char *) malloc(sizeof(char) * count);
    }
    MPI_Gatherv((void *) output.c_str(), outputSize, MPI_CHAR, combinedOutput, counts, displacements, MPI_CHAR, 0,
                MPI_COMM_WORLD);
/*
    if (rank == 0 && fsave) {
        fprintf(fsave, combinedOutput);
    }
    */
}

void initSystem() {
    int sendCount[maxRank];
    int sendDisplacement[maxRank];

    // Allocate space for particles
    particle_t *particles = (particle_t *) malloc(globalParticleCount * sizeof(particle_t));
    particle_t *particlesToSend = (particle_t *) malloc(globalParticleCount * sizeof(particle_t));
    ownedParticles = (particle_t *) malloc(globalParticleCount * sizeof(particle_t));

    // Initialize grid and world (for all processors)
    grid_init(sqrt(0.0005 * globalParticleCount) / 0.01);
    set_size(globalParticleCount);
    cellsPerProcess = gridColumns * (int) ceil(gridColumns / (double) maxRank);

    printf("CellsPerProcess: %d\n", cellsPerProcess);

    // Initialize and prepare particles at the master node
    if (rank == 0) {
        initializeAtRoot(particles, particlesToSend, sendCount, sendDisplacement);
    }

    // Distribute the particle count to all processors
    MPI_Scatter(sendCount, 1, MPI_INT, &maxPosition, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Distribute the actual particles
    MPI_Scatterv(particlesToSend, sendCount, sendDisplacement, particleType, ownedParticles, maxPosition,
                 particleType,
                 0, MPI_COMM_WORLD);

    // Initialize world zones.
    // TODO We allocate quite a bit of memory here. Does the density provide us with any guarantees?
    borrowedLower.particles = (particle_t *) malloc(sizeof(particle_t) * globalParticleCount);
    borrowedUpper.particles = (particle_t *) malloc(sizeof(particle_t) * globalParticleCount);
    borrowedLower.coordinateStart = (cellsPerProcess * rank) - gridColumns;
    borrowedUpper.coordinateStart = (cellsPerProcess * (rank + 1));

    ownedLower.particles = (particle_t *) malloc(sizeof(particle_t) * globalParticleCount);
    ownedUpper.particles = (particle_t *) malloc(sizeof(particle_t) * globalParticleCount);
    ownedLower.coordinateStart = cellsPerProcess * rank;
    ownedUpper.coordinateStart = min((cellsPerProcess * (rank + 1)) - gridColumns, gridColumns * (gridColumns - 1));

    // Initialize grid with received particles
    for (int i = 0; i < maxPosition; ++i) {
        grid_add(&ownedParticles[i]);
    }

#ifdef DEBUG
    WHEN_DEBUGGING(validateScatter(particlesToSend, ownedParticles, maxPosition));
    assert(cellsPerProcess % gridColumns == 0);
    assert(borrowedLower.coordinateStart % gridColumns == 0);
    assert(borrowedUpper.coordinateStart % gridColumns == 0);
    assert(ownedLower.coordinateStart % gridColumns == 0);
    assert(ownedUpper.coordinateStart % gridColumns == 0);

    if (rank == 0) {
        printf("---- World zones ----\n");
        printf("Borrowed lower: %d\n", borrowedLower.coordinateStart);
        printf("Owned lower: %d\n", ownedLower.coordinateStart);
        printf("Owned upper: %d\n", ownedUpper.coordinateStart);
        printf("Borrowed upper: %d\n", borrowedUpper.coordinateStart);
        printf("---- World zones ----\n");
    }
#endif

    free(particles);
    free(particlesToSend);
}

void validateScatter(particle_t *particlesToSend, particle_t *localParticles, int localCount) {
    if (rank == 0) {
        for (int i = 0; i < localCount; i++) {
            particle_t &copied = particlesToSend[i];
            particle_t &particle = localParticles[i];

            assert(copied.ax == particle.ax && copied.ay == particle.ay && copied.vx == particle.vx
                   && copied.vy == particle.vy && copied.x == particle.x && copied.y == particle.y);
        }
        printf("Local variables match\n");
    }
}

void initializeAtRoot(particle_t *particles, particle_t *particlesToSend, int *sendCount, int *sendDisplacement) {
    int sendStart = 0;
    int currentProcessor = 0;
    int counter = 0;

    init_particles(globalParticleCount, particles);

    for (int i = 0; i < globalParticleCount; ++i) {
        grid_add(&particles[i]);
    }
    for (int i = 0; i < gridColumns * gridColumns; i++) {
        auto cell = grid_get_at(i);
        for (auto particle : cell) {
            memcpy(particlesToSend + counter, particle, sizeof(particle_t));
            counter++;
        }
        if ((i > 0 && i % cellsPerProcess == 0) || i == (gridColumns * gridColumns) - 1) {
            sendCount[currentProcessor] = counter - sendStart;
            sendDisplacement[currentProcessor] = sendStart;
            sendStart = counter;
            currentProcessor++;
        }
    }
    grid_reset();

#ifdef DEBUG
    validateRootInitialization(currentProcessor, counter, particlesToSend);
#endif
}

void validateRootInitialization(int currentProcessor, int counter, particle_t *particlesToSend) {
    assert(currentProcessor == maxRank);
    assert(counter == globalParticleCount);

    counter = 0;
    for (int i = 0; i < gridColumns * gridColumns; i++) {
        auto cell = grid_get_at(i);
        for (auto particle : cell) {
            particle_t &copied = particlesToSend[counter];
            assert(copied.ax == particle->ax && copied.ay == particle->ay && copied.vx == particle->vx
                   && copied.vy == particle->vy && copied.x == particle->x && copied.y == particle->y);
            counter++;
        }
    }
}

void initMPI(int &argc, char **&argv) {//  set up MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &maxRank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Type_contiguous(6, MPI_DOUBLE, &particleType);
    MPI_Type_commit(&particleType);
}

int printHelp() {
    printf("Options:\n");
    printf("-h to see this help\n");
    printf("-n <int> to set the number of particles\n");
    printf("-o <filename> to specify the output file name\n");
    return 0;
}

void wait_for_debugger() {
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i) {
        sleep(5);
    }
}