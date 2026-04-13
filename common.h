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
static int read_particles_from_input_file(const char *input_file, Particle *particles, int num_particles) {
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

static int read_particles_from_output_file(const char *input_file, Particle *particles, int num_particles) {
    FILE *fp = fopen(input_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error opening input file: %s\n", input_file);
        return 1;
    }

    // read file + populate particles
    for (int i = 0; i < num_particles; i++) {
        float x, y;
        int q;
        float p;
        if (fscanf(fp, "%f,%f,%d,%f", &x, &y, &q, &p) != 4) {
            fprintf(stderr, "Error reading particle data from file: %s\n", input_file);
            fclose(fp);
            return 1;
        }
        particles[i].z = x + y*I;
        particles[i].x = x;
        particles[i].y = y;
        particles[i].q = q;
        particles[i].p = p;
    }
    fclose(fp);
    return 0;
}

static int compare_particles(const void *a, const void *b) {
    Particle *p1 = (Particle *)a;
    Particle *p2 = (Particle *)b;

    if (p1->x < p2->x) return -1;
    if (p1->x > p2->x) return 1;
    if (p1->y < p2->y) return -1;
    if (p1->y > p2->y) return 1;
    return 0; // particles are at the same position
}

static int write_particles_to_output_file(const char *output_file, Particle *particles, int num_particles) {
    // sort particles by x-coordinate then y-coordinate to ensure easy reading and comparison across files
    qsort(particles, num_particles, sizeof(Particle), compare_particles);

    FILE *fp = fopen(output_file, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error opening output file: %s\n", output_file);
        return 1;
    }

    for (int i = 0; i < num_particles; i++) {
        fprintf(fp, "%.8f,%.8f,%d,%.8f\n", particles[i].x, particles[i].y, particles[i].q, particles[i].p);
    }

    fclose(fp);
    return 0;
}

static int write_particles_to_data_file(const char *output_file, Particle *particles, int num_particles) {
    FILE *fp = fopen(output_file, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error opening output file: %s\n", output_file);
        return 1;
    }

    for (int i = 0; i < num_particles; i++) {
        fprintf(fp, "%.8f,%.8f,%d\n", particles[i].x, particles[i].y, particles[i].q);
    }

    fclose(fp);
    return 0;
}