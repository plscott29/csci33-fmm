#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"

/* intended usage:
    ./compare <output_file1> <output_file2> <num_particles>
*/
int main(int argc, char * argv[]) {
    // read inputs from command line
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <output_file1> <output_file2> <num_particles>\n", argv[0]);
        fprintf(stderr, "output_file1: path to the first output file containing particle data\n");
        fprintf(stderr, "output_file2: path to the second output file containing particle data\n");
        fprintf(stderr, "num_particles: number of particles\n");
        return 1;
    }

    char *output_file1 = argv[1];
    char *output_file2 = argv[2];
    int num_particles = atoi(argv[3]);

    // read particle data from the first output file
    Particle *particles_fmm = (Particle *)malloc(num_particles * sizeof(Particle));
    if (particles_fmm == NULL) {
        fprintf(stderr, "Error allocating memory for particles\n");
        return 1;
    }

    if (read_particles_from_output_file(output_file1, particles_fmm, num_particles) != 0) {
        fprintf(stderr, "Error reading particles from file\n");
        free(particles_fmm);
        return 1;
    }

    // read particle data from the second output file to compare against the first output file
    Particle *particles_direct = (Particle *)malloc(num_particles * sizeof(Particle));
    if (particles_direct == NULL) {
        fprintf(stderr, "Error allocating memory for direct particles\n");
        free(particles_fmm);
        return 1;
    }
    if (read_particles_from_output_file(output_file2, particles_direct, num_particles) != 0) {
        fprintf(stderr, "Error reading particles from file\n");
        free(particles_fmm);
        free(particles_direct);
        return 1;
    }

    // compare potentials from FMM and direct calculations normalized by  A=max(|p|) across all particles
    float A = 0.0f;
    for (int i = 0; i < num_particles; i++) {
        A += fabsf(particles_fmm[i].q);  // charge is the same across FMM and direct
    }

    float max_error = 0.0f;
    for (int i = 0; i < num_particles; i++) {
        float error = fabsf(particles_fmm[i].p - particles_direct[i].p);
        if (error > max_error) {
            max_error = error;
        }
    }
    
    float epsilon = max_error / A;
    printf("Relative precision epslon between %s and %s: %.5e\n", output_file1, output_file2, epsilon);

    // free allocated memory
    free(particles_fmm);
    free(particles_direct);
    return 0;
}
