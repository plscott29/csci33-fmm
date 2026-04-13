#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"

#define CHARGE_MIN -5
#define CHARGE_MAX 5

/* intended usage:
    ./data_writer <num_particles> <output_file>
*/
int main(int argc, char * argv[]) {
    // read inputs from command line
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_particles> <output_file>\n", argv[0]);
        fprintf(stderr, "num_particles: number of particles\n");
        fprintf(stderr, "output_file: path to the output file where results will be written\n");
        return 1;
    }

    int num_particles = atoi(argv[1]);
    char *output_file = argv[2];

    // allocate particles array and initialize with random data
    Particle *particles = (Particle *)malloc(num_particles * sizeof(Particle));
    if (particles == NULL) {
        fprintf(stderr, "Error allocating memory for particles\n");
        return 1;
    }

    // seed random number generator for reproducibility
    srand(time(NULL));

    for (int i = 0; i < num_particles; i++) {
        particles[i].x = (float)rand() / ((float)RAND_MAX + 1.0f);
        particles[i].y = (float)rand() / ((float)RAND_MAX + 1.0f);
        particles[i].q = CHARGE_MIN + rand() % (CHARGE_MAX - CHARGE_MIN + 1);
    }

    // write particles to output file
    if (write_particles_to_data_file(output_file, particles, num_particles) != 0) {
        fprintf(stderr, "Error writing particles to file\n");
        free(particles);
        return 1;
    }

    free(particles);
    return 0;
}