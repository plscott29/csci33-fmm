#include <complex.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"

int compute_direct_potentials_parallel(Particle *particles, int num_particles, int num_threads) {
    // TODO: parallelize direct potential calculation using OpenMP, Newton's 3rd Law, and memory tiling
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

    if (read_particles_from_file(input_file, particles, num_particles) != 0) {
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
    if (write_particles_to_file(output_file, particles, num_particles) != 0) {
        fprintf(stderr, "Error writing particles to file\n");
        free(particles);
        return 1;
    }

    // free allocated memory
    free(particles);
    return 0;
}
