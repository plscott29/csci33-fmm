#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

typedef struct {
    float x;
    float y;
    float complex z;
    int q;
    float p;
} Particle;

typedef struct {
    //int idx;
    int level;
    int start;
    int end;
    float x_mid;
    float y_mid;
    float complex z_mid;
    float complex expansions[]; // to store local & multipole expansions, size determined by 2p
} Node;

typedef struct {
    int start;
    int first_quad;
    int second_quad;
    int third_quad;
    int end;
} Partition;

// outputs the index of the first child given a parent index
// the other children are the following three indices
int child_idx(int parent_idx) {return (parent_idx+1)*4;}
//outputs the index of the parent given a child index
int parent_idx(int child_idx) {return child_idx/4-1;}

//float mid_x(int node_idx) {}

//sort particles, return bounds for the quadrants as a Partition
Partition four_sort(Particle *particles, Node *node) { //int start, int end, float x_mid, float y_mid) {
    Particle temp_particle;
    int low_idx = node->start;
    int high_idx = node->end-1;
    // sort on x values first
    for (int i =0; i < node->end - node->start; i++) {
        if (particles[low_idx].x <= node->x_mid)
            low_idx++;
        else {
            if (particles[high_idx].x > node->x_mid) {
                high_idx--;
            }
            else { //swap low and high idx particles
                temp_particle = particles[high_idx];
                particles[high_idx] = particles[low_idx];
                particles[low_idx] = temp_particle;
                low_idx++;
            }
        }
    }

    // now [0, low_idx] <= x_mid and [low_idx,end-start] > x_mid
    int second_quad = low_idx;

    // sort on y values for all points with x <= x_mid
    low_idx = node->start;
    high_idx = second_quad-1;
    for (int i =0; i < second_quad - node->start; i++) {
        if (particles[low_idx].y <= node->y_mid)
            low_idx++;
        else {
            if (particles[high_idx].y > node->y_mid) {
                high_idx--;
            }
            else {
                temp_particle = particles[high_idx];
                particles[high_idx] = particles[low_idx];
                particles[low_idx] = temp_particle;
                low_idx++;
            }
        }
    }

    // now [0,low_idx] <= y_mid and [low_idx, second_quad] > y_mid
    int first_quad = low_idx;

    // sort on y values for all points with x > x_mid
    low_idx = second_quad;
    high_idx = node->end - 1;
    for (int i =0; i < node->end - second_quad; i++) {
        if (particles[low_idx].y <= node->y_mid)
            low_idx++;
        else {
            if (particles[high_idx].y > node->y_mid) {
                high_idx--;
            }
            else {
                temp_particle = particles[high_idx];
                particles[high_idx] = particles[low_idx];
                particles[low_idx] = temp_particle;
                low_idx++;
            }
        }
    }

    // now [0,low_idx] is leq y_mid and [low_idx, second_quad] is gt y_mid
    int third_quad = low_idx;

    Partition quadrants = {node->start, first_quad, second_quad, third_quad, node->end};
    return quadrants;
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
    Node *root = malloc(sizeof(Node));
    //root->idx = -1;
    root->level = 0;
    root->start = 0; root->end = num_particles;
    root->x_mid = 0.5; root->y_mid = 0.5;
    root->z_mid = 0.5 + 0.5*I;
    //printf("root->z_mid r %f, i%f\n", creal(root->z_mid), cimag(root->z_mid));
    // need to manually populate the first four elements of nodes[] using four_sort
    Partition boundaries = four_sort(particles, root);

    // manually set the first level
    nodes[0].level = 1; nodes[0].start = boundaries.start; nodes[0].end = boundaries.first_quad;
    nodes[0].x_mid = 0.25; nodes[0].y_mid = 0.25; nodes[0].z_mid = nodes[0].x_mid + nodes[0].y_mid*I;
    nodes[1].level = 1; nodes[1].start = boundaries.first_quad; nodes[1].end = boundaries.second_quad;
    nodes[1].x_mid = 0.25; nodes[1].y_mid = 0.75; nodes[1].z_mid = nodes[1].x_mid + nodes[1].y_mid*I;
    nodes[2].level = 1; nodes[2].start = boundaries.second_quad; nodes[2].end = boundaries.third_quad;
    nodes[2].x_mid = 0.75; nodes[2].y_mid = 0.25; nodes[2].z_mid = nodes[2].x_mid + nodes[2].y_mid*I;
    nodes[3].level = 1; nodes[3].start = boundaries.third_quad; nodes[3].end = boundaries.end;
    nodes[3].x_mid = 0.75; nodes[3].y_mid = 0.75; nodes[3].z_mid = nodes[3].x_mid + nodes[3].y_mid*I;

    // set the rest of the levels
    int p_idx = 0;
    int c_idx;
    int num_boxes = 4;
    for (int level = 1; level < num_levels; level++) { // level iteration number refers to the level of the parent
        for (int i = 0; i < num_boxes; i++) { // num_boxes is the number of parent boxes
            // sort and get the partition of the subarray of particles
            boundaries = four_sort(particles, &nodes[p_idx]);

            // set four nodes
            //printf("quads: %d, %d, %d, %d, %d\n", boundaries.start, boundaries.first_quad,
            //    boundaries.second_quad, boundaries.third_quad, boundaries.end);
            c_idx = child_idx(p_idx);
            float shift = 1.0/pow(2, level+2);

            nodes[c_idx].level = level+1;
            nodes[c_idx].start = boundaries.start; nodes[c_idx].end = boundaries.first_quad;
            nodes[c_idx].x_mid = nodes[p_idx].x_mid - shift; nodes[c_idx].y_mid = nodes[p_idx].y_mid - shift;
            nodes[c_idx].z_mid = nodes[c_idx].x_mid + nodes[c_idx].y_mid*I;
            c_idx++;

            nodes[c_idx].level = level+1;
            nodes[c_idx].start = boundaries.first_quad; nodes[c_idx].end = boundaries.second_quad;
            nodes[c_idx].x_mid = nodes[p_idx].x_mid - shift; nodes[c_idx].y_mid = nodes[p_idx].y_mid + shift;
            nodes[c_idx].z_mid = nodes[c_idx].x_mid + nodes[c_idx].y_mid*I;
            c_idx++;

            nodes[c_idx].level = level+1;
            nodes[c_idx].start = boundaries.second_quad; nodes[c_idx].end = boundaries.third_quad;
            nodes[c_idx].x_mid = nodes[p_idx].x_mid + shift; nodes[c_idx].y_mid = nodes[p_idx].y_mid - shift;
            nodes[c_idx].z_mid = nodes[c_idx].x_mid + nodes[c_idx].y_mid*I;
            c_idx++;

            nodes[c_idx].level = level+1;
            nodes[c_idx].start = boundaries.third_quad; nodes[c_idx].end = boundaries.end;
            nodes[c_idx].x_mid = nodes[p_idx].x_mid + shift; nodes[c_idx].y_mid = nodes[p_idx].y_mid + shift;
            nodes[c_idx].z_mid = nodes[c_idx].x_mid + nodes[c_idx].y_mid*I;
            //c_idx++; not used

            p_idx++;
        }

        //start_idx = (start_idx+1)*4;
        num_boxes = num_boxes*4;
    }
    //printf("node %d after calling construct_tree: r %f, i %f\n", 82, creal(nodes[82].expansions[1]), cimag(nodes[82].expansions[1]));

    free(root);
    return 0;
}

void ofs(Particle *particles, Node *node, int P) {
    for (int  j = node->start; j < node->end; j++) {
        //printf("Charge number = %d\n", j);
        node->expansions[0] = node->expansions[0] + particles[j].q*1.0;
        float complex offset = particles[j].z - node->z_mid;
        //printf("particles[j].z r %f, i %f\n", creal(particles[j].z), cimag(particles[j].z));
        //printf("node->z_mid r %f, i %f\n", creal(node->z_mid), cimag(node->z_mid));
        //printf("offset r %f, i %f\n", creal(offset), cimag(offset));
        for (int k=1; k < P; k++) {
            //printf("Multipole expansion k = %d, real %f imag %f\n",
            //    k, creal(cpow(offset, k)*particles[j].q/k), cimag(cpow(offset, k)*particles[j].q/k));
            node->expansions[k] = node->expansions[k] - cpow(offset, k)*particles[j].q/k;
        }
    }
}

int calculate_multipole(Particle *particles, Node *nodes, int num_particles, int num_levels, int P) {
    int num_nodes = 0; // total nodes
    for (int k = 1; k <= num_levels; k++) {num_nodes = num_nodes + pow(4,k);}
    int num_leaves = pow(4, num_levels);
    for (int i = num_nodes-1; num_nodes-num_leaves-1 <= i; i--) {
        //printf("node number: %d\n", i);
        //printf("node %d before: r %f, i %f\n", i, creal(nodes[i].expansions[1]), cimag(nodes[i].expansions[1]));
        //printf("next node (%d) before: r %f, i %f\n", i-1, creal(nodes[i-1].expansions[1]), cimag(nodes[i-1].expansions[1]));
        ofs(particles, &nodes[i], P);
        //printf("node %d after: r %f, i %f\n\n", i, creal(nodes[i].expansions[1]), cimag(nodes[i].expansions[1]));
    }
    return 0;
}

/* intended usage:
    ./fmm <input_file> <output_file> <num_levels> <p> <num_particles> <num_threads>
*/
int main(int argc, char * argv[])
{
    // read inputs from command line
    if (argc != 7) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> <num_levels> <p> <num_particles> <num_threads>\n", argv[0]);
        fprintf(stderr, "input_file: path to the input file containing particle data\n");
        fprintf(stderr, "output_file: path to the output file where results will be written\n");
        fprintf(stderr, "num_levels: depth of FMM expansion tree\n");
        fprintf(stderr, "P: order of multipole expansion\n");
        fprintf(stderr, "num_particles: number of particles\n");
        fprintf(stderr, "num_threads: number of threads to use\n");
        return 1;
    }

    char *input_file = argv[1];
    char *output_file = argv[2];
    int num_levels = atoi(argv[3]);
    int P = atoi(argv[4]);
    int num_particles = atoi(argv[5]);
    int num_threads = atoi(argv[6]);

    // read input particle data from file
    Particle *particles = (Particle *)malloc(num_particles * sizeof(Particle));
    //Particle *particles = (Particle *)calloc(num_particles, sizeof(Particle));
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
        particles[i].z = particles[i].x + particles[i].y*I; //TODO: incorporate this into the fscanf line above
    }
    fclose(fp);

    // allocate memory for tree structure
    // num_nodes = \sum_{l=0}^{num_levels} 4^l
    int num_nodes = 0; // total nodes
    for (int k = 1; k <= num_levels; k++) {num_nodes = num_nodes + pow(4,k);}

    Node *nodes = (Node *)malloc(num_nodes * (sizeof(Node) + 2*P*sizeof(float complex)));
    if (nodes == NULL) {
        fprintf(stderr, "Error allocating memory for tree nodes\n");
        free(particles);
        return 1;
    }

    // run step 1: tree construction & sorting
    if (construct_tree(particles, nodes, num_particles, num_levels) != 0) {
        fprintf(stderr, "Error constructing tree\n");
        free(particles);
        free(nodes);
        return 1;
    }

    // run step 2: calculate multipole expansions
    if (calculate_multipole(particles, nodes, num_particles, num_levels, P) != 0) {
        fprintf(stderr, "Error calculating multipole expansions\n");
        free(particles);
        free(nodes);
        return 1;
    }

    printf("Node idx: %d, mpole expansion 1: %f, %f\n", 83, creal(nodes[83].expansions[1]), cimag(nodes[83].expansions[1]));

    return 0;
}