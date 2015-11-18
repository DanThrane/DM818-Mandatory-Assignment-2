#include <mpi.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sstream>

#include "grid.h"
#include "mpi.h"

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
int localCount;

/**
 * The particles owned by this processor
 */
particle_t *ownedParticles;

GhostZone ownedUpper;
GhostZone ownedLower;

GhostZone borrowedUpper;
GhostZone borrowedLower;

/**
 * Grid coordinate for when the borrowed from below ghost zone starts
 */
int borrowedGhostZoneLowerStart;

int main(int argc, char **argv) {
    printf("This is the new one\n");
    //  process command line parameters
    globalParticleCount = read_int(argc, argv, "-n", 1000);
    char *savename = read_string(argc, argv, "-o", NULL);
    FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : NULL;
    if (find_option(argc, argv, "-h") >= 0) { return printHelp(); }

    initMPI(argc, argv);
    initSystem();

    //  simulate a number of time steps
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {
        MPI_Barrier(MPI_COMM_WORLD);
        // TODO Merge
        // Purge the ghost zone from last iteration
        grid_purge(borrowedGhostZoneUpperStart, borrowedGhostZoneUpperStart + gridColumns);
        // Send and receive new ghost zones (TODO Do in parallel)
        if (rank > 0) {
            MPI_Recv(borrowedGhostZoneUpper, prevGhostCount, particleType, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < maxRank - 1) {
            MPI_Send(ownedGhostZone, ownedGhostZoneParticleCount, particleType, rank + 1, 0, MPI_COMM_WORLD);
        }
        // Add the new ghost zone to the grid
        if (rank > 0) {
            for (int i = 0; i < prevGhostCount; i++) {
                grid_add(&borrowedGhostZoneUpper[i]);
            }
        }

        //  compute forces
        for (int i = 0; i < localCount; i++) {
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

        //  move particles
        //  and update their position in the grid
        for (int i = 0; i < localCount; i++) {
            grid_remove(&ownedParticles[i]);
            move(ownedParticles[i]);
            grid_add(&ownedParticles[i]);
        }

        //  save if necessary
        if (savename && step % SAVEFREQ == 0) {
            for (int i = 0; i < localCount; i++) {
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

    if (rank == 0) { printf("'done'\n"); }

    // Release resources
    // TODO Release stuff
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
        printf("Mallocing %d characters\n", count);
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
    MPI_Scatter(sendCount, 1, MPI_INT, &localCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("I'm %d and I will receive %d\n", rank, localCount);
    // Distribute the actual particles
    MPI_Scatterv(particlesToSend, sendCount, sendDisplacement, particleType, ownedParticles, localCount, particleType,
                 0, MPI_COMM_WORLD);

    // Initialize world zones.
    // TODO We allocate quite a bit of memory here. Does the density provide us with any guarantees?
    borrowedLower.particles = (particle_t *) malloc(sizeof(particle_t) * globalParticleCount);
    borrowedUpper.particles = (particle_t *) malloc(sizeof(particle_t) * globalParticleCount);
    borrowedLower.coordinateStart = (cellsPerProcess * rank) - gridColumns;
    borrowedUpper.coordinateStart = (cellsPerProcess * rank + 1);

    ownedLower.particles = (particle_t *) malloc(sizeof(particle_t) * globalParticleCount);
    ownedUpper.particles = (particle_t *) malloc(sizeof(particle_t) * globalParticleCount);
    ownedLower.coordinateStart = cellsPerProcess * rank;
    ownedUpper.coordinateStart = (cellsPerProcess * rank + 1) - gridColumns;

    // Initialize grid with received particles
    for (int i = 0; i < localCount; ++i) {
        grid_add(&ownedParticles[i]);
    }

#ifdef DEBUG
    validateScatter(particlesToSend, ownedParticles, localCount);
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