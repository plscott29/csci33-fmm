#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

/* Define structs: Particle, Node, Partition */
typedef struct {
    float complex z;
    float x;
    float y;
    int q;
    float p;
} Particle;

typedef struct {
    int level;
    int start;
    int end;
    float x_mid;
    float y_mid;
    float complex z_mid;
    float complex expansions[]; // to store local & multipole expansions, size determined by 2p
} Node;

typedef struct {
    int quadrant_bounds[4][2];  // Quadrant ordering: SW/NW/SE/NE
} Partition;


/* Helper functions:
    - child_idx: given a parent node index, returns the index of the first child in the quadtree
    - parent_idx: given a child node index, returns the index of the parent node in the quadtree
    - choose: computes the binomial coefficient "n choose k"
*/
int child_idx(int parent_idx) { return (parent_idx+1)*4; }
int parent_idx(int child_idx) { return child_idx/4-1; }

int choose(int n, int k) {
    int c = 1;
    
    // nCk = nC(n-k) --> exploit symmetry & use smaller of (k, n-k) for computation
    if (k > n-k) { k=n-k; }

    // compute nCk = n! / (k! * (n-k)!)
    for (int i = n; i > k; i--)     { c = c*i; }    // n*(n-1)*...*(k+1)
    for (int i = 1; i <= n-k; i++)  { c = c/i; }    // divide by (n-k)!
    
    return c;
}

/* 
    input:
        - particles: array of particles
        - node: node to partition, which points to a subarray of particles with start & end indices

    output:
        - returns a Partition struct with the indices of four quadrants of input node (inclusive)

    Sorts subarray of particles corresponding to the input node into four quadrants based on x and y values
*/
Partition four_sort(Particle *particles, Node *node) {
    Particle tmp;
    int low = node->start;
    int high = node->end;

    // sort particles by x values first
    while (low < high) {
        if (particles[low].x <= node->x_mid)
            low++;
        else if (particles[high].x > node->x_mid)
            high--;
        else {
            // swap particles at indices low & high 
            tmp = particles[high];
            particles[high] = particles[low];
            particles[low] = tmp;
            low++;
            high--;
        }
    }

    // particles with x <= x_mid are now in range [node->start, low-1]
    // and particles with x > x_mid are in range [low, node->end]
    int x_split = low;  // index of first particle in right half (x > x_mid)

    // sort particles by y values for those with x <= x_mid
    // (i.e., all particles belonging to left half)
    low = node->start;
    high = x_split - 1;
    while (low < high) {
        if (particles[low].y <= node->y_mid)
            low++;
        else if (particles[high].y > node->y_mid)
            high--;
        else {
            tmp = particles[high];
            particles[high] = particles[low];
            particles[low] = tmp;
            low++;
            high--;
        }
    }

    // particles with x <= x_mid and y <= y_mid are now in range [node->start, low-1]
    // and particles with x <= x_mid and y > y_mid are in range [low, x_split-1]
    int left_y_split = low;

    // sort particles by y values for those with x > x_mid
    low = x_split;
    high = node->end;
    while (low < high) {
        if (particles[low].y <= node->y_mid)
            low++;
        else if (particles[high].y > node->y_mid)
            high--;
        else {
            tmp = particles[high];
            particles[high] = particles[low];
            particles[low] = tmp;
            low++;
            high--;
        }
    }

    // particles with x > x_mid and y <= y_mid are now in range [x_split, low-1]
    // and particles with x > x_mid and y > y_mid are in range [low, node->end]
    int right_y_split = low;

    return (Partition) {
        .quadrant_bounds = {
            {node->start, left_y_split-1},
            {left_y_split, x_split-1},
            {x_split, right_y_split-1},
            {right_y_split, node->end}
        }
    };
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
    // initialize root node
    Node *root = malloc(sizeof(Node));
    if (root == NULL) {
        return -1;
    }

    // root->idx = -1; 
    root->level = 0;
    root->start = 0; root->end = num_particles-1;
    
    // TODO: parameterize bounding box?
    root->x_mid = 0.5; root->y_mid = 0.5; 
    root->z_mid = 0.5 + 0.5*I;

    // sort particles into four quadrants at the root
    Partition partitions = four_sort(particles, root);
    
    // quadrant ordering: SW, NW, SE, NE
    float x_mids[4] = {0.25, 0.25, 0.75, 0.75};
    float y_mids[4] = {0.25, 0.75, 0.25, 0.75};

    // create first level of child nodes
    for (int i=0; i < 4; i++) {
        nodes[i].level = 1;
        nodes[i].start = partitions.quadrant_bounds[i][0];
        nodes[i].end = partitions.quadrant_bounds[i][1];
        nodes[i].x_mid = x_mids[i];
        nodes[i].y_mid = y_mids[i];
        nodes[i].z_mid = nodes[i].x_mid + nodes[i].y_mid*I;
    }

    // set the rest of the levels
    int p_idx = 0;
    int c_idx;
    int num_boxes = 4;

    for (int level = 1; level < num_levels; level++) {  // 'level' refers to level of parent
        for (int i = 0; i < num_boxes; i++) {           // 'num_boxes' is number of parent boxes
            // sort particles and get the partitions of subarray of particles for this parent node
            partitions = four_sort(particles, &nodes[p_idx]);

            // set four nodes corresponding to the four quadrants of this parent node
            c_idx = child_idx(p_idx);

            // calculate offset from parent midpoint to child midpoint:
            // shift = (1/2^(level+1)) / 2 = 1/2^(level+2)
            float shift = 1.0 / pow(2, level+2);

            // quadrant offsets from parent to child with ordering: SW, NW, SE, NE
            float x_mid_shifts[4] = {-shift, -shift, shift, shift};
            float y_mid_shifts[4] = {-shift, shift, -shift, shift};

            for (int q = 0; q < 4; q++) {
                nodes[c_idx + q].level = level+1;
                nodes[c_idx + q].start = partitions.quadrant_bounds[q][0];
                nodes[c_idx + q].end = partitions.quadrant_bounds[q][1];
                nodes[c_idx + q].x_mid = nodes[p_idx].x_mid + x_mid_shifts[q];
                nodes[c_idx + q].y_mid = nodes[p_idx].y_mid + y_mid_shifts[q];
                nodes[c_idx + q].z_mid = nodes[c_idx + q].x_mid + nodes[c_idx + q].y_mid*I;
            }

            p_idx++;
        }

        num_boxes = num_boxes*4;
    }

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

void ofo(Node *nodes, int p_idx, int P, int *binom) {
    int c_idx = child_idx(p_idx);
    for (int i = 0; i < 4; i++) {
        float complex offset = nodes[c_idx+i].z_mid - nodes[p_idx].z_mid;
        for (int r = 1; r <= P; r++) {
            float complex offset_power = 1;
            for (int s = r; s >= 0; s--) { //iterate backwards to fill in multipole expansions from low degree to high degree
                nodes[p_idx].expansions[r-1] += binom[(r-1)*P+(s-1)]*offset_power*nodes[c_idx+i].expansions[s-1];
                offset_power *= offset;
            }
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

    //precalculate binomial coefs
    int *binom = (int*)malloc(P*P*sizeof(int));
    for (int i = 1; i <= P; i++) {
        for (int j=1;j <= i; j++) {
            binom[P*(i-1)+(j-1)] = choose(i,j);
        }
    }

    for (int i = num_nodes-num_leaves-2; i > 0; i--) {
        printf("node number: %d\n", i);
        printf("node %d before: r %f, i %f\n", i, creal(nodes[i].expansions[1]), cimag(nodes[i].expansions[1]));
        printf("next node (%d) before: r %f, i %f\n", i-1, creal(nodes[i-1].expansions[1]), cimag(nodes[i-1].expansions[1]));
        ofo(nodes, i, P, binom);
        printf("node %d after: r %f, i %f\n\n", i, creal(nodes[i].expansions[1]), cimag(nodes[i].expansions[1]));
    }

    free(binom);
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
    if (particles == NULL) {
        fprintf(stderr, "Error allocating memory for particles\n");
        return 1;
    }

    FILE *fp = fopen(input_file, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error opening input file: %s\n", input_file);
        free(particles);
        return 1;
    }

    for (int i = 0; i < num_particles; i++) {
        float x, y;
        int q;
        if (fscanf(fp, "%f,%f,%d", &x, &y, &q) != 3) {
            fprintf(stderr, "Error reading particle data from file: %s\n", input_file);
            free(particles);
            fclose(fp);
            return 1;
        }
        particles[i].z = x + y*I;
        particles[i].x = x;
        particles[i].y = y;
        particles[i].q = q;
    }
    fclose(fp);

    // allocate memory for tree structure: \sum_{l=0}^{num_levels} 4^l
    int num_nodes = 0;
    for (int k = 1; k <= num_levels; k++) { num_nodes = num_nodes + pow(4,k); }

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