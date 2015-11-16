#include <mpi.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#include "grid.h"
#include "common.h"

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

    MPI_Datatype PARTICLE;
    MPI_Type_contiguous(6, MPI_DOUBLE, &PARTICLE);
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

    int prevGhostCount;
    particle_t *prevGhostZone;

    int nextGhostCount;
    particle_t *nextGhostZone;


#ifdef DEBUG
    if (rank == 0) {
        for (int i = 0; i < localCount; i++) {
            particle_t &copied = particlesToSend[i];
            particle_t &particle = localParticles[i];

            assert(copied.ax == particle.ax && copied.ay == particle.ay && copied.vx == particle.vx
                   && copied.vy == particle.vy && copied.x == particle.x && copied.y == particle.y);
        }
        printf("Local variables match\n");
    }
#endif

    for (int i = 0; i < localCount; ++i) {
        grid_add(&localParticles[i]);
    }

/*

    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {
        for (int i = 0; i < n; i++) {

            particles[i].ax = particles[i].ay = 0;

            // traverse included neighbors
            for (int offsetX = -1; offsetX <= 1; offsetX++) {
                for (int offsetY = -1; offsetY <= 1; offsetY++) {
                    const std::vector<particle_t *> &cell =
                            grid_get_collisions_at_neighbor(&particles[i], offsetX, offsetY);

                    for (auto particle : cell) {
                        apply_force(particles[i], *particle);
                    }
                }
            }
        }

        //  move particles
        //  and update their position in the grid
        for (int i = 0; i < n; i++) {
            grid_remove(&particles[i]);
            move(particles[i]);
            grid_add(&particles[i]);
        }

        //  save if necessary
        if (fsave && (step % SAVEFREQ) == 0) {
            save(fsave, n, particles);
        }
    }
    simulation_time = read_timer() - simulation_time;

    if (rank == 0) {
        printf("n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time);
    }
*/
    if (rank == 0) { printf("'done'\n"); }
    // Release resources
    // TODO Release stuff
    if (fsave) {
        fclose(fsave);
    }

    MPI_Finalize();
    return 0;
}
