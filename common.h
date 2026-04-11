#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

/* define Particle struct that will be used across FMM and direct potential calculations */
typedef struct {
    float complex z;
    float x;
    float y;
    int q;
    float p;
} Particle;

/* define helper functions used across FMM and direct potential calculations */
static int read_particles_from_file(const char *input_file, Particle *particles, int num_particles) {
    FILE *fp = fopen(input_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error opening input file: %s\n", input_file);
        return 1;
    }

    // read file + populate particles
    for (int i = 0; i < num_particles; i++) {
        float x, y;
        int q;
        if (fscanf(fp, "%f,%f,%d", &x, &y, &q) != 3) {
            fprintf(stderr, "Error reading particle data from file: %s\n", input_file);
            fclose(fp);
            return 1;
        }
        particles[i].z = x + y*I;
        particles[i].x = x;
        particles[i].y = y;
        particles[i].q = q;
    }
    fclose(fp);
    return 0;
}

static int write_particles_to_file(const char *output_file, Particle *particles, int num_particles) {
    FILE *fp = fopen(output_file, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error opening output file: %s\n", output_file);
        return 1;
    }
    for (int i = 0; i < num_particles; i++) {
        fprintf(fp, "%f,%f,%d,%f\n", particles[i].x, particles[i].y, particles[i].q, particles[i].p);
    }
    fclose(fp);
    return 0;
}