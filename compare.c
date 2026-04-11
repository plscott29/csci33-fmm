#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"

/* intended usage:
    ./compare <fmm_output_file> <direct_output_file> <num_particles>
*/
int main(int argc, char * argv[]) {
    // read inputs from command line
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <fmm_output_file> <direct_output_file> <num_particles>\n", argv[0]);
        fprintf(stderr, "fmm_output_file: path to the FMM output file containing particle data\n");
        fprintf(stderr, "direct_output_file: path to the direct output file where results will be written\n");
        fprintf(stderr, "num_particles: number of particles\n");
        return 1;
    }

    char *fmm_output_file = argv[1];
    char *direct_output_file = argv[2];
    int num_particles = atoi(argv[3]);

    // read particle data from FMM output file
    Particle *particles_fmm = (Particle *)malloc(num_particles * sizeof(Particle));
    if (particles_fmm == NULL) {
        fprintf(stderr, "Error allocating memory for particles\n");
        return 1;
    }

    if (read_particles_from_file(fmm_output_file, particles_fmm, num_particles) != 0) {
        fprintf(stderr, "Error reading particles from file\n");
        free(particles_fmm);
        return 1;
    }

    // read particle data from direct output file to compare against FMM results
    Particle *particles_direct = (Particle *)malloc(num_particles * sizeof(Particle));
    if (particles_direct == NULL) {
        fprintf(stderr, "Error allocating memory for direct particles\n");
        free(particles_fmm);
        return 1;
    }
    if (read_particles_from_file(direct_output_file, particles_direct, num_particles) != 0) {
        fprintf(stderr, "Error reading direct particles from file\n");
        free(particles_fmm);
        free(particles_direct);
        return 1;
    }

    // compare potentials from FMM and direct calculations
    float max_error = 0.0f;
    for (int i = 0; i < num_particles; i++) {
        float error = fabsf(particles_fmm[i].p - particles_direct[i].p);
        if (error > max_error) {
            max_error = error;
        }
    }
    printf("Maximum error between FMM and brute-force potentials: %e\n", max_error);

    // free allocated memory
    free(particles_fmm);
    free(particles_direct);
    return 0;
}
