#include <complex.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"

#define TILE_SIZE 64

int compute_direct_potentials_parallel(Particle *particles, int num_particles, int num_threads) {
    /* parallelize direct potential calculation using OpenMP, Newton's 3rd Law, and memory tiling */
    omp_set_num_threads(num_threads);

    int num_tiles = ceil((float) num_particles / TILE_SIZE);

    // initialize potentials to 0 for each particle
    #pragma omp parallel for
    for (int i = 0; i < num_particles; i++) {
        particles[i].p = 0.0f;
    }

    // diagonal tiles (independent and fully parallelizable)
    #pragma omp parallel for
    for (int tile = 0; tile < num_tiles; tile++) {
        int start = tile * TILE_SIZE;
        int end = fmin(start + TILE_SIZE, num_particles);

        for (int p_idx = start; p_idx < end; p_idx++) {
            Particle *p = &particles[p_idx];
            for (int q_idx = p_idx + 1; q_idx < end; q_idx++) {
                Particle *q = &particles[q_idx];
                float complex offset = p->z - q->z;
                float complex log_offset = clogf(offset);
                float contribution_p = crealf(q->q * log_offset);
                float contribution_q = p->q * crealf(log_offset);
                
                p->p += contribution_p;
                q->p += contribution_q;
            }
        }
    }

    // off-diagonal tiles (parallelize over tiles, but not within tiles)
    for (int r = 0; r < num_tiles - 1; r++) {
        #pragma omp parallel for
        for (int s = 0; s < num_tiles / 2; s++) {
            int ti, tj;

            if (s == 0) { // slot 0: (T-1, r)
                ti = num_tiles - 1;
                tj = r;
            } else { // slot 1: (r+1, T-2), generally: ((r+s) % (T-1), (r-s+T-1) % (T-1))
                ti = (r + s) % (num_tiles-1);
                tj = (r + num_tiles - s) % (num_tiles-1);
            }

            // process tile pair (ti, tj) with Newton's 3rd Law
            int start_i = ti * TILE_SIZE;
            int end_i = fmin(start_i + TILE_SIZE, num_particles);
            int start_j = tj * TILE_SIZE;
            int end_j = fmin(start_j + TILE_SIZE, num_particles);
            for (int p_idx = start_i; p_idx < end_i; p_idx++) {
                Particle *p = &particles[p_idx];
                for (int q_idx = start_j; q_idx < end_j; q_idx++) {
                    Particle *q = &particles[q_idx];
                    float complex offset = p->z - q->z;
                    float complex log_offset = clogf(offset);
                    float contribution_p = crealf(q->q * log_offset);
                    float contribution_q = p->q * crealf(log_offset);
                    
                    p->p += contribution_p;
                    q->p += contribution_q;
                }
            }
        } // implicit barrier: all threads finish round r before starting round r+1 (no tile is accessed by two rounds simultaneously)
    }

    return 0;
}

int compute_direct_potentials(Particle *particles, int num_particles) {
    for (int p_idx = 0; p_idx < num_particles; p_idx++) {
        Particle *p = &particles[p_idx];
        p->p = 0.0f;
        for (int q_idx = 0; q_idx < num_particles; q_idx++) {
            if (p_idx != q_idx) {
                Particle q = particles[q_idx];
                float complex offset = p->z - q.z;
                p->p += crealf(q.q * clogf(offset));
            }
        }
    }
    return 0;
}

/* intended usage:
    ./direct <input_file> <output_file> <num_particles> <is_parallel> <num_threads>
*/
int main(int argc, char * argv[]) {
    // read inputs from command line
    if (argc < 6) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> <num_particles> <is_parallel> <num_threads>\n", argv[0]);
        fprintf(stderr, "input_file: path to the input file containing particle data\n");
        fprintf(stderr, "output_file: path to the output file where results will be written\n");
        fprintf(stderr, "num_particles: number of particles\n");
        fprintf(stderr, "is_parallel: whether to run in parallel (1) or sequential (0)\n");
        fprintf(stderr, "num_threads: number of threads to use\n");
        return 1;
    }

    char *input_file = argv[1];
    char *output_file = argv[2];
    int num_particles = atoi(argv[3]);
    int is_parallel = atoi(argv[4]);
    int num_threads = atoi(argv[5]);

    // read input particle data from file
    Particle *particles = (Particle *)malloc(num_particles * sizeof(Particle));
    if (particles == NULL) {
        fprintf(stderr, "Error allocating memory for particles\n");
        return 1;
    }

    if (read_particles_from_input_file(input_file, particles, num_particles) != 0) {
        fprintf(stderr, "Error reading particles from file\n");
        free(particles);
        return 1;
    }

    // direct particle-wise potentials on particles array
    if (is_parallel) {
        compute_direct_potentials_parallel(particles, num_particles, num_threads);
    } else {
        compute_direct_potentials(particles, num_particles);
    }

    // write results to output file
    if (write_particles_to_output_file(output_file, particles, num_particles) != 0) {
        fprintf(stderr, "Error writing particles to file\n");
        free(particles);
        return 1;
    }

    // free allocated memory
    free(particles);
    return 0;
}
