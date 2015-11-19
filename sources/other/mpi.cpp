#include <mpi.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sstream>

#include "grid.h"
#include "mpi.h"
#include "common.h"

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

int lowerBorder;
int upperBorder;

GhostZone ownedUpper;
GhostZone ownedLower;

GhostZone borrowedUpper;
GhostZone borrowedLower;

std::vector<particle_t *> insertionsUpper;
std::vector<particle_t *> insertionsLower;

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
    printf("Exchanging particles between %d and %d\n", rank, rank + (1 * multiplier));
    MPI_Sendrecv(&owned.particleCount, 1, MPI_INT, rank + (1 * multiplier), 0, &borrowed.particleCount, 1,
                 MPI_INT, rank + (1 * multiplier), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    printf("[%d] I will receive %d from %d\n", rank, borrowed.particleCount, rank + (1 * multiplier));
    printf("[%d] I will send %d to %d\n", rank, owned.particleCount, rank + (1 * multiplier));

    // Exchange particles
    MPI_Sendrecv(owned.particles, owned.particleCount, particleType, rank + (1 * multiplier), 0,
                 borrowed.particles, borrowed.particleCount, particleType, rank + (1 * multiplier), 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    VALIDATE_GHOST_ZONE(owned);
    VALIDATE_GHOST_ZONE(borrowed);

    printf("Finished exchanging particles between %d and %d\n", rank, rank + (1 * multiplier));
}

particle_t *prepareInsertions(std::vector<particle_t *> insertions) {
    int i = 0;
    particle_t *result = (particle_t *) malloc(sizeof(particle_t) * insertions.size());
    for (auto particle : insertions) {
        memcpy(&result[i], particle, sizeof(particle_t));
        i++;
    }
    return result;
}

particle_t *exchangeInsertions(std::vector<particle_t *> &insertions, int multiplier, int *outputCount) {
    int count;
    int sendingCount = (int) insertions.size();

    particle_t *prepared = prepareInsertions(insertions);
    particle_t *receiveBuffer;

    MPI_Sendrecv(&sendingCount, 1, MPI_INT, rank + (1 * multiplier), 0, &count, 1, MPI_INT, rank + (1 * multiplier), 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    receiveBuffer = (particle_t *) malloc(sizeof(particle_t) * count);
    printf("[%d] Receiving %d inserted particles!\n", rank, count);
    MPI_Sendrecv(prepared, sendingCount, particleType, rank + (1 * multiplier), 0, receiveBuffer, count, particleType,
                 rank + (1 * multiplier), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    *outputCount = count;
    free(prepared);
    return receiveBuffer;
}

void updateGrid(GhostZone &zone, std::vector<particle_t *> &localInsertions) {
    VALIDATE_GHOST_ZONE(zone);

    grid_purge(zone.coordinateStart, zone.coordinateStart + gridColumns);
    for (int i = 0; i < zone.particleCount; i++) {
        grid_add(&zone.particles[i]);
    }
    for (auto particle : localInsertions) {
        grid_add(particle);
    }

    VALIDATE_GHOST_ZONE(zone);
}

void mergeInsertedInLocallyOwned(int insertedCount, particle_t *inserted) {
    int start = maxPosition;

    // The incoming buffer is a temporary one. Move the particles to the owned buffer
    memcpy(&ownedParticles[maxPosition], inserted, sizeof(particle_t) * insertedCount);
    maxPosition += insertedCount;

    for (int i = start; i < start + insertedCount; i++) {
        grid_add(&ownedParticles[i]);
    }
}

void exchangeInformationAbove(particle_t **insertedLower, int *insertedLowerCount) {
    if (rank < maxRank - 1) {
        exchangeParticles(ownedUpper, borrowedUpper, 1);
        *insertedLower = exchangeInsertions(insertionsUpper, 1, insertedLowerCount);
    }
}

void exchangeInformationBelow(particle_t **insertedUpper, int *insertedUpperCount) {
    if (rank > 0) {
        exchangeParticles(ownedLower, borrowedLower, -1);
        *insertedUpper = exchangeInsertions(insertionsLower, -1, insertedUpperCount);
    }
}

void exchangeInformationWithNeighborHood() {
    int insertedUpperCount = 0;
    particle_t *insertedUpper = NULL;
    int insertedLowerCount = 0;
    particle_t *insertedLower = NULL;

    // Prepare ghost zones such that we're ready to exchange them with our neighbors
    prepareGhostZoneForExchange(ownedUpper);
    prepareGhostZoneForExchange(ownedLower);

    if (rank % 2 == 0) {
        // Even ranks communicate with processors above us first
        exchangeInformationAbove(&insertedLower, &insertedLowerCount);
        exchangeInformationBelow(&insertedUpper, &insertedUpperCount);
    } else {
        // Even ranks communicate with processors below us first
        exchangeInformationBelow(&insertedUpper, &insertedUpperCount);
        exchangeInformationAbove(&insertedLower, &insertedLowerCount);
    }

    // Disable tracking while we update our ghost zones
    grid_disable_track();

    // First insert the complete state changes from our neighbor merged with out own local insertions in these zones
    updateGrid(borrowedUpper, insertionsUpper);
    updateGrid(borrowedLower, insertionsLower);

    // Then we get the insertions that occurred into our owned zones from our neighbors
    mergeInsertedInLocallyOwned(insertedLowerCount, insertedLower);
    mergeInsertedInLocallyOwned(insertedUpperCount, insertedUpper);

    grid_enable_track();

    // Clear insertions from last iteration
    insertionsLower.clear();
    insertionsUpper.clear();

    // Release inserted
    if (insertedLower == NULL) free(insertedLower);
    if (insertedUpper == NULL) free(insertedUpper);
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

    grid_track_insertions(borrowedLower.coordinateStart, borrowedLower.coordinateStart + gridColumns, &insertionsLower);
    grid_track_insertions(borrowedUpper.coordinateStart, borrowedUpper.coordinateStart + gridColumns, &insertionsUpper);

    //  simulate a number of time steps
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {
        if (rank == 0) printf("Step %d\n", step);
        MPI_Barrier(MPI_COMM_WORLD);
        WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));
        exchangeInformationWithNeighborHood();
        WHEN_DEBUGGING(validate_grid(0, __FILE__, __LINE__));

        //  compute forces
        for (int i = 0; i < maxPosition; i++) {
            ownedParticles[i].ax = ownedParticles[i].ay = 0;

            // traverse included neighbors
            for (int offsetX = -1; offsetX <= 1; offsetX++) {
                for (int offsetY = -1; offsetY <= 1; offsetY++) {
                    const std::vector<particle_t *> &cell =
                            grid_get_collisions_at_neighbor(&ownedParticles[i], offsetX, offsetY);

                    for (auto particle : cell) {
                        apply_force(ownedParticles[i], *particle);
                    }
                }
            }
        }

        // Move and update particles
        for (int i = 0; i < maxPosition; i++) {
            particle_t *particle = &ownedParticles[i];
            grid_remove(particle);
            move(*particle);
            grid_add(particle);

            int coordinate = get_particle_coordinate(particle);
            if (coordinate >= lowerBorder && coordinate < upperBorder) {
                // Out of bounds - Particle will be purged soon. Take last element and move it here
                maxPosition--;
                if (i != maxPosition) {
                    memcpy(&ownedParticles[i], &ownedParticles[maxPosition], sizeof(particle_t));
                }
            }
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

    if (rank == 0 && fsave) {
        fprintf(fsave, combinedOutput);
    }
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
    cellsPerProcess = (gridColumns * gridColumns) / maxRank;

    // Initialize and prepare particles at the master node
    if (rank == 0) {
        initializeAtRoot(particles, particlesToSend, sendCount, sendDisplacement);
    }

    // Distribute the particle count to all processors
    MPI_Scatter(sendCount, 1, MPI_INT, &maxPosition, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Distribute the actual particles
    MPI_Scatterv(particlesToSend, sendCount, sendDisplacement, particleType, ownedParticles, maxPosition, particleType,
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
    ownedUpper.coordinateStart = (cellsPerProcess * (rank + 1)) - gridColumns;

    lowerBorder = ownedLower.coordinateStart;
    upperBorder = ownedUpper.coordinateStart + gridColumns;

    // Initialize grid with received particles
    for (int i = 0; i < maxPosition; ++i) {
        grid_add(&ownedParticles[i]);
    }

#ifdef DEBUG
    validateScatter(particlesToSend, ownedParticles, maxPosition);
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