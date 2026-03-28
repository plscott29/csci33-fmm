#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct {
    float x;
    float y;
    int q;
    float p;
} Particle;

typedef struct {
    int idx;
    int level;
    int start;
    int end;
    float expansions[]; // to store local & multipole expansions, size determined by 2p
} Node;

/* intended usage:
    ./fmm <input_file> <output_file> <num_levels> <p> <num_particles> <num_threads>
*/
int main(int argc, char * argv[])
{
    fprintf(stdout, "Hello, World!\n");

    // read inputs from command line
    if (argc != 7) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> <num_levels> <p> <num_particles> <num_threads>\n", argv[0]);
        fprintf(stderr, "input_file: path to the input file containing particle data\n");
        fprintf(stderr, "output_file: path to the output file where results will be written\n");
        fprintf(stderr, "num_levels: depth of FMM expansion tree\n");
        fprintf(stderr, "p: order of multipole expansion\n");
        fprintf(stderr, "num_particles: number of particles\n");
        fprintf(stderr, "num_threads: number of threads to use\n");
        return 1;
    }

    char *input_file = argv[1];
    char *output_file = argv[2];
    int num_levels = atoi(argv[3]);
    int p = atoi(argv[4]);
    int num_particles = atoi(argv[5]);
    int num_threads = atoi(argv[6]);

    // read input particle data from file
    Particle *particles = (Particle *)malloc(num_particles * sizeof(Particle));
    FILE *fp = fopen(input_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error opening input file: %s\n", input_file);
        free(particles);
        return 1;
    }
    for (int i = 0; i < num_particles; i++) {
        if (fscanf(fp, "%f,%f,%d", &particles[i].x, &particles[i].y, &particles[i].q) != 3) {
            fprintf(stderr, "Error reading particle data from file: %s\n", input_file);
            free(particles);
            fclose(fp);
            return 1;
        }
    }
    fclose(fp);

    // allocate memory for tree structure
    // num_nodes = \sum_{l=0}^{num_levels} 4^l
    int num_nodes = ((1 << (2 * num_levels)) - 1) / 3; // total nodes
    Node *nodes = (Node *)malloc(num_nodes * (sizeof(Node) + 2*p*sizeof(float)));
    if (nodes == NULL) {
        fprintf(stderr, "Error allocating memory for tree nodes\n");
        free(particles);
        return 1;
    }

    // run step 1: tree construction & sorting
    /*
    if (construct_tree(particles, nodes, num_particles, num_levels) != 0) {
        fprintf(stderr, "Error constructing tree\n");
        free(particles);
        free(nodes);
        return 1;
    }
    */
    return 0;
}

/*
    input: 
        - particles: array of particles
        - nodes: pre-allocated array of tree nodes
        - num_particles: number of particles
        - num_levels: depth of the tree
    output:
        - returns 0 on success, non-zero on failure

    In-place sort input array of particles and populate every node
    in the tree with (i) idx, (ii) level, (iii) start, and (iv) end.
*/
int construct_tree(Particle *particles, Node *nodes, int num_particles, int num_levels)
{
    return 0;
}