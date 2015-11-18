#include <mpi.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <sstream>

#include "grid.h"
#include "common.h"

std::ostringstream localOutput;

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

//  benchmarking program
int main(int argc, char **argv) {
    //  process command line parameters
    if (find_option(argc, argv, "-h") >= 0) {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);
    char *savename = read_string(argc, argv, "-o", NULL);

    //  set up MPI
    int n_proc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //  allocate generic resources
    FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : NULL;
    particle_t *particles = (particle_t *) malloc(n * sizeof(particle_t));
    particle_t *particlesToSend = (particle_t *) malloc(n * sizeof(particle_t));
    particle_t *localParticles = (particle_t *) malloc(n * sizeof(particle_t));

    printf("fsave %s\n", savename);

    MPI_Datatype PARTICLE;
    MPI_Type_contiguous(6, MPI_DOUBLE, &PARTICLE);
    MPI_Type_contiguous(1, MPI_INT, &PARTICLE);
    MPI_Type_commit(&PARTICLE);

    int sendCount[n_proc];
    int sendDisplacement[n_proc];

    grid_init(sqrt(0.0005 * n) / 0.01);
    int cellsPerProcess = (gridColumns * gridColumns) / n_proc;
    int currentProcessor = 0;
    int sendStart = 0;

    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    set_size(n);
    if (rank == 0) {
//        wait_for_debugger();
        init_particles(n, particles);

        for (int i = 0; i < n; ++i) {
            grid_add(&particles[i]);
        }

        int counter = 0;
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
        assert(currentProcessor == n_proc);
        assert(counter == n);

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
        printf("done copying shit\n");
        for (int i = 0; i < n_proc; i++) {
            printf("Process %d should receive %d\n", i, sendCount[i]);
        }
#endif
    }

    int localCount;
    MPI_Scatter(sendCount, 1, MPI_INT, &localCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("I am %d and I will receive %d\n", rank, localCount);
    MPI_Scatterv(particlesToSend, sendCount, sendDisplacement, PARTICLE, localParticles, localCount, PARTICLE, 0,
                 MPI_COMM_WORLD);

    int prevGhostCount = gridColumns;
    particle_t *prevGhostZone = (particle_t *) malloc(sizeof(particle_t) * prevGhostCount);
    int prevGhostStart = (cellsPerProcess * rank) - gridColumns;
    int prevGhostEnd = cellsPerProcess * rank;

    int nextGhostCount = gridColumns;
    particle_t *ownedGhostZone = &localParticles[localCount - gridColumns];


#ifdef DEBUG
    if (rank == 0) {
        for (int i = 0; i < localCount; i++) {
            particle_t &copied = particlesToSend[i];
            particle_t &particle = localParticles[i];

            assert(copied.ax == particle.ax && copied.ay == particle.ay && copied.vx == particle.vx
                   && copied.vy == particle.vy && copied.x == particle.x && copied.y == particle.y
                    && copied.id == particle.id);
        }
        printf("Local variables match\n");
    }
#endif

    for (int i = 0; i < localCount; ++i) {
        grid_add(&localParticles[i]);
    }

    //  simulate a number of time steps
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {
        MPI_Barrier(MPI_COMM_WORLD);
        // Purge the ghost zone from last iteration
        grid_purge(prevGhostStart, prevGhostEnd);
        // Send and receive new ghost zones (TODO Do in parallel)
        if (rank > 0) {
            MPI_Recv(prevGhostZone, prevGhostCount, PARTICLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < n_proc - 1) {
            MPI_Send(ownedGhostZone, nextGhostCount, PARTICLE, rank + 1, 0, MPI_COMM_WORLD);
        }
        // Add the new ghost zone to the grid
        if (rank > 0) {
            for (int i = 0; i < prevGhostCount; i++) {
                grid_add(&prevGhostZone[i]);
            }
        }

        //  compute forces
        for (int i = 0; i < localCount; i++) {
            localParticles[i].ax = localParticles[i].ay = 0;

            // traverse included neighbors
            for (int offsetX = -1; offsetX <= 1; offsetX++) {
                for (int offsetY = -1; offsetY <= 1; offsetY++) {
                    const std::vector<particle_t *> &cell =
                            grid_get_collisions_at_neighbor(&localParticles[i], offsetX, offsetY);

                    for (auto particle : cell) {
                        apply_force(localParticles[i], *particle);
                    }
                }
            }

        }

        //  move particles
        //  and update their position in the grid
        for (int i = 0; i < localCount; i++) {
            grid_remove(&localParticles[i]);
            move(localParticles[i]);
            grid_add(&localParticles[i]);
        }

        //  save if necessary
        if (savename && step % SAVEFREQ == 0) {
            for (int i = 0; i < localCount; i++) {
                localOutput << localParticles[i].x << " " << localParticles[i].y << "\n";
            }
        }
        // Communicate to rank 0 with result
    }
    simulation_time = read_timer() - simulation_time;

    if (savename) {
        int counts[n_proc];
        char *combinedOutput;
        int *displacements;

        int count = 0;
        std::string output = localOutput.str();
        int outputSize = (int) output.size();

        MPI_Gather(&outputSize, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            displacements = (int *) malloc(sizeof(int) * n_proc);
            for (int i = 0; i < n_proc; i++) {
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

    if (rank == 0) { printf("n = %d, simulation time = %g seconds\n", n, simulation_time); }

    if (rank == 0) { printf("'done'\n"); }
    // Release resources
    // TODO Release stuff
    if (fsave) {
        fclose(fsave);
    }

    MPI_Finalize();
    return 0;
}
