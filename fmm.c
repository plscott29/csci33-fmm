#include <complex.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

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
    float complex *expansions; // to store local & multipole expansions, size determined by 2p
} Node;

typedef struct {
    int quadrant_bounds[4][2];  // Quadrant ordering: SW/NW/SE/NE
} Partition;

typedef struct {
    int nodes[3][3];
} Neighborhood;

typedef struct {
    int nodes[27];
} WellSeparated;


/* Helper functions:
    - child_idx: given a parent node index, returns the index of the first child in the quadtree
    - parent_idx: given a child node index, returns the index of the parent node in the quadtree
    - choose: computes the binomial coefficient "n choose k"
*/
int child_idx(int parent_idx) { return (parent_idx * 4) + 1; }
int parent_idx(int child_idx) { return (child_idx - 1) / 4; }
int level_start_idx(int level) {
    int idx = 0;
    for (int k = 0; k < level; k++) { idx = idx + pow(4,k); }
    return idx;
}

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
    Node *root = &nodes[0];                 // TODO: optimize by not storing root?
    root->level = 0;
    root->start = 0; root->end = num_particles-1;
    root->x_mid = 0.5; root->y_mid = 0.5;   // TODO: parameterize bounding box?
    root->z_mid = 0.5 + 0.5*I;

    // set the rest of the levels
    int p_idx = 0;
    int c_idx;
    int num_boxes = 1;     // 'num_boxes' is number of parent boxes

    for (int level = 0; level < num_levels; level++) {  // 'level' refers to level of parent node(s)
        for (int i = 0; i < num_boxes; i++) {
            // sort particles and get the partitions of subarray of particles for this parent node
            Partition partitions = four_sort(particles, &nodes[p_idx]);

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

    return 0;
}

/*
    input:
        - particles: array of Particle structs containing the particle positions and charges
        - node: pointer to the Node for which the multipole expansion is being computed
        - P: number of terms in the multipole expansion
    output:
        - node->expansions: array of complex numbers of length P in which the computed
          multipole expansion coefficients for this node are stored

    Compute "outgoing from source" expansions for all leaf nodes in the tree.
    For each particle in the node's particle range [node->start, node->end], compute
    the contribution of that particle to the multipole expansion around the node's center
    'z_mid' and accumulate it into node->expansions according to:
        M_0 = M_0 + q_j
        M_p = M_p - q_j/p * (z_j - node->z_mid)^p for p = 1, ..., P-1
    where z_j is the particle position and q_j is the particle charge.
*/
void ofs(Particle *particles, Node *node, int P) {
    for (int j = node->start; j <= node->end; j++) {
        node->expansions[0] = node->expansions[0] + (float) particles[j].q;

        float complex offset = particles[j].z - node->z_mid;
        for (int p = 1; p < P; p++) {
            node->expansions[p] = node->expansions[p] - cpow(offset, p) * (float) particles[j].q / (float) p;
        }
    }
}

/*
    input:
        - nodes: tree of Nodes representing FMM hiearchical decomposition
        - node_idx: index of node in nodes array for which multipole expansion is being computed
        - P: number of terms in multipole expansion
        - binom: precomputed array of binomial coefficients (of size P*P)
    output:
        - nodes[node_idx].expansions: multipole expansion coefficients for this node are updated
          in place using multipole expansions of its children.
          
    Compute "outgoing from outgoing" multipole expansions for this node according to:
        M_r = sum_{s=0}^{r} binom(r,s) * (z_c - z_p)^(r-s) * M_s^c
    where M_r are the multipole expansion coefficients of the parent node,
    z_c is the center of the child node c, z_p is the center of the parent node,
    r is the index of the multipole expansion coefficient of the parent node,
    s is the index of the multipole expansion coefficient of the child node c,
    M_s^c are the multipole expansion coefficients of the child node c,
    and binom(r,s) is the binomial coefficient "r choose s".
*/
void ofo(Node *nodes, int node_idx, int P, int *binom) {
    int c_idx = child_idx(node_idx);    // index of first child of this node

    // iterate over children
    for (int i = 0; i < 4; i++) {
        float complex offset = nodes[c_idx+i].z_mid - nodes[node_idx].z_mid;

        for (int r = 1; r <= P; r++) {
            // accumulate powers of offset for use in binomial expansion
            float complex offset_power_acc = 1;

            // iterate backwards to fill in multipole expansions from low degree to high degree
            for (int s = r; s > 0; s--) {
                // accumulate contribution of child multipole expansion to parent
                int r_choose_s = binom[(r-1)*P + (s-1)];
                nodes[node_idx].expansions[r-1] += r_choose_s * offset_power_acc * nodes[c_idx+i].expansions[s-1];
                offset_power_acc *= offset;
            }
        }
    }
}

/*
    input:
        - nodes: tree of Nodes representing FMM hiearchical decomposition
        - node_idx: index of node in nodes array for which incoming-from-incoming
          multipole expansion is being computed
        - P: number of terms in multipole expansion
        - binom: precomputed array of binomial coefficients (of size P*P)
    output:
        - nodes[node_idx].expansions[P..2*P-1]: incoming expansion coefficients for this node
          are updated in place using incoming expansions from the parent node

    Compute "incoming from incoming" expansion for this node by translating the incoming
    expansions from the parent node according to:
        u_\sigma += T^ifi_\sigma,\tau u_\tau
    where u_\sigma are the incoming expansion coefficients for this node, u_\tau are the incoming
    expansion coefficients of the parent node, and T^ifi_\sigma,\tau is the upper-triangular
    translation operator that maps incoming expansions from the parent to the child node.

*/
void ifi(Node *nodes, int node_idx, int P, int *binom) {
    float complex *incoming_expansions = &nodes[node_idx].expansions[P];

    int p_idx = parent_idx(node_idx);
    float complex offset =  nodes[node_idx].z_mid - nodes[p_idx].z_mid; // child - parent

    // iterate from highest degree of incoming expansion down to lowest degree
    for (int r = P; r > 0; r--) {
        // accumulate powers of offset for use in binomial expansion
        float complex offset_power_acc = 1;

        // iterate backwards to fill in multipole expansions from low degree to high degree
        for (int p = r; p <= P; p++) {
            // accumulate contribution of child multipole expansion to parent
            int p_choose_r = binom[(p-1)*P + (r-1)];
            incoming_expansions[r-1] += p_choose_r * offset_power_acc * nodes[p_idx].expansions[P+(p-1)];
            offset_power_acc *= offset;
        }
    }
}

int search_for_neighbor(Node *nodes, int node_idx, char direction, float distance) {
    int idx_lower_bound = level_start_idx(nodes[node_idx].level);
    int idx_upper_bound = level_start_idx(nodes[node_idx].level+1)-1;

    float test_y_dist;
    float test_x_dist;
    bool aligned = 0;
    float ratio;
    int iteration = 0;
    int jump_idx;
    int test_idx = node_idx;

    switch (direction) {
        case 'n':
            while (test_idx <= idx_upper_bound) {
                jump_idx = pow(4, iteration)-level_start_idx(iteration); // +4-1=3, +16-4-1=11, +64-16-4-1=43, etc.
                test_idx = node_idx + jump_idx; // we want to see if this node is the north neighbor
                printf("node %d, jump %d, test %d\n", node_idx, jump_idx, test_idx);
                test_y_dist = cimag(nodes[test_idx].z_mid - nodes[node_idx].z_mid); printf("test distance: %f\n", test_y_dist);
                ratio = test_y_dist/distance; printf("ratio: %f\n", ratio);

                aligned = creal(nodes[test_idx].z_mid) == creal(nodes[node_idx].z_mid); printf("aligned: %d\n", aligned);

                if (0.8 < ratio && ratio < 1.2 && aligned) {
                    return test_idx;
                }
                iteration++;
            }
            break;

        case 's':
            while (test_idx >= idx_lower_bound) {
                jump_idx = pow(4, iteration)-level_start_idx(iteration);
                test_idx = node_idx - jump_idx;
                test_y_dist = cimag(nodes[node_idx].z_mid - nodes[test_idx].z_mid);
                ratio = test_y_dist/distance;
                aligned = creal(nodes[test_idx].z_mid) == creal(nodes[node_idx].z_mid);
                if (0.95 < ratio && ratio < 1.05 && aligned) {
                    return test_idx;
                }
                iteration++;
            }
            break;

        case 'e':
            while (test_idx <= idx_upper_bound) {
                jump_idx = 2*(pow(4, iteration)-level_start_idx(iteration));
                test_idx = node_idx + jump_idx;
                test_x_dist = creal(nodes[test_idx].z_mid - nodes[node_idx].z_mid);
                ratio = test_x_dist/distance;
                aligned = cimag(nodes[test_idx].z_mid) == cimag(nodes[node_idx].z_mid);
                if (0.8 < ratio && ratio < 1.2 && aligned) {
                    return test_idx;
                }
                iteration++;
            }
            break;

        case 'w':
            while (test_idx >= idx_lower_bound) {
                jump_idx = 2*(pow(4, iteration)-level_start_idx(iteration));
                test_idx = node_idx - jump_idx;
                //printf("node %d, jump %d, test %d\n", node_idx, jump_idx, test_idx);
                test_x_dist = creal(nodes[node_idx].z_mid - nodes[test_idx].z_mid);
                ratio = test_x_dist/distance;
                aligned = cimag(nodes[test_idx].z_mid) == cimag(nodes[node_idx].z_mid);
                if (0.95 < ratio && ratio < 1.05 && aligned) {
                    return test_idx;
                }
                iteration++;
            }
            break;
    }

    return -1;
}

/*
    input:
        - nodes: tree of Nodes representing FMM hiearchical decomposition
        - node_idx: index of node in nodes array for which incoming-from-outgoing
          multipole expansion is being computed
        - P: number of terms in multipole expansion
        - binom: precomputed array of binomial coefficients (of size P*P)
    output:
        - nodes[node_idx].expansions[P..2*P-1]: incoming expansion coefficients for this node
          are updated in place using outgoing expansions from nodes in the interaction list
    
    Compute "incoming from outgoing" expansion for this node from nodes in the interation list
    according to:
        u_\sigma += \sum_{c in interaction list} T^ifo_\sigma,c q_c
    where u_\sigma are the incoming expansions for this node, q_c are outgoing
    expansion coefficients of nodes in the interaction list of this node, and T^ifo_\sigma,c
    is the translation operator that maps outgoing expansions of the interaction list nodes
    to incoming expansions at this node
*/
void ifo(Node *nodes, int node_idx, int P) {
    float complex *incoming_expansions = &nodes[node_idx].expansions[P];

    int interaction_list[] = {0};       // TODO: construct interaction list for each node
    for (int i = 0; i < -1; i++) {      // TODO: currently set so iteration doesn't occur
        int c_idx = interaction_list[i];
        float complex offset = nodes[node_idx].z_mid - nodes[c_idx].z_mid;          // current - interaction list node
        float complex reverse_offset = nodes[c_idx].z_mid - nodes[node_idx].z_mid;  // interaction list node - current

        // r = 0 case
        incoming_expansions[0] += nodes[c_idx].expansions[0] * clog(offset);
        float complex offset_power_acc = -1*reverse_offset; // (c_sigma - c_tau)^p starting at p=1
        for (int p = 1; p < P; p++) {
            incoming_expansions[0] += nodes[c_idx].expansions[p] / offset_power_acc;
            offset_power_acc *= -1*reverse_offset; // (c_sigma - c_tau)^p
        }

        // r > 0 case
        offset_power_acc = reverse_offset; // (c_sigma - c_tau)^r starting at r=1
        for (int r = 1; r < P; r++) {
            incoming_expansions[r] -= nodes[c_idx].expansions[0] / (r * offset_power_acc);

            float complex offset_power_acc_inner = offset_power_acc; // (c_sigma - c_tau)^r+p starting at (c_sigma - c_tau)^r
            for (int p = 1; p < P; p++) {
                offset_power_acc_inner *= -1*reverse_offset; // (c_sigma - c_tau)^(r+p)

                int binom_factor = choose(r+p-1, p-1);  // (r + p - 1) choose (p - 1)
                incoming_expansions[r] += nodes[c_idx].expansions[p] * binom_factor / offset_power_acc_inner;
            }

            offset_power_acc *= reverse_offset; // (c_sigma - c_tau)^r
        }
    }
}

int calculate_multipole(Particle *particles, Node *nodes, int num_particles, int num_levels, int num_nodes, int P) {
    /* 
     * Upwards pass: compute outgoing expansion of each box from leaf nodes up to the root. 
     * For a leaf node, compute expansion directly particles in the node ("outgoing from sources"); 
     * for a parent node, use outgoing expansions of its children to ("outgoing from outgoing")
     */

    // compute multipole expansions for all leaf nodes: "outgoing from sources"
    int num_leaves = pow(4, num_levels);
    for (int i = 0; i < num_leaves; i++) {
        int leaf_idx = num_nodes - num_leaves + i;
        ofs(particles, &nodes[leaf_idx], P);
    }

    // pre-compute binomial co-efficients
    int *binom = (int*)malloc(P*P*sizeof(int));
    for (int r = 1; r <= P; r++) {
        for (int s = 1; s <= P; s++) {
            binom[P*(r-1)+(s-1)] = choose(r,s);
        }
    }

    // compute multipole expansions for all non-leaf nodes: "outgoing from outgoing"
    // iterate over all non-leaf nodes in reverse order (from bottom of tree up to root)
    for (int i = (num_nodes - num_leaves) - 1; i > 0; i--) {
        ofo(nodes, i, P, binom);
    }

    free(binom);
    return 0;
}

int calculate_local(Particle *particles, Node *nodes, int num_particles, int num_levels, int num_nodes, int P) {
    /*
     * Downwards pass: compute incoming expansion for every node in a pass over all nodes,
     * going from larger->smaller. For each node, combine the incoming expansion of its parent
     * with the contributions from the outgoing expansions of all boxes in its interaction list.
    */

    // set incoming expansion of root & level 1 nodes to zero
    for (int i = 0; i < (1 + 4); i++) { // root node + 4 level 1 nodes
        for (int p = 0; p < P; p++) {
            memset(&nodes[i].expansions[P], 0, P*sizeof(float complex));
        }
    }

    // pre-compute binomial co-efficients
    int *binom = (int*)malloc(P*P*sizeof(int));
    for (int r = 1; r <= P; r++) {
        for (int s = 1; s <= P; s++) {
            binom[P*(r-1)+(s-1)] = choose(r,s);
        }
    }

    // iterate from level l=2 to leaf nodes, setting incoming expansions
    for (int l = 2; l <= num_levels; l++) {
        int nodes_in_level = pow(4, l);
        int level_start_idx = (pow(4, l) - 1) / 3;
        for (int i = 0; i < nodes_in_level; i++) {
            int node_idx = level_start_idx + i;
            ifi(nodes, node_idx, P, binom);
            ifo(nodes, node_idx, P);
        }
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

    // read file + populate particles
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
    int num_nodes = level_start_idx(num_levels+1);

    // allocate space for nodes
    Node *nodes = (Node *)malloc(num_nodes * (sizeof(Node)));
    if (nodes == NULL) {
        fprintf(stderr, "Error allocating memory for tree nodes\n");
        free(particles);
        return 1;
    }
    for (int i=0; i < num_nodes; i++) {
        nodes[i].expansions = (float complex *) calloc(2 * P, sizeof(float complex));
        if (nodes[i].expansions == NULL) {
            fprintf(stderr, "Error allocating memory for expansions of node %d\n", i);
            for (int j = 0; j < i; j++) {
                free(nodes[j].expansions);
            }
            free(nodes);
            free(particles);
            return 1;
        }
    }

    // run step 1: tree construction & sorting
    if (construct_tree(particles, nodes, num_particles, num_levels) != 0) {
        fprintf(stderr, "Error constructing tree\n");
        free(particles);
        for (int i = 0; i < num_nodes; i++) {
            free(nodes[i].expansions);
        }
        free(nodes);
        return 1;
    }

    // run step 2: calculate multipole expansions (upwards pass)
    if (calculate_multipole(particles, nodes, num_particles, num_levels, num_nodes, P) != 0) {
        fprintf(stderr, "Error calculating multipole expansions\n");
        free(particles);
        for (int i = 0; i < num_nodes; i++) {
            free(nodes[i].expansions);
        }
        free(nodes);
        return 1;
    }

    // run step 3: calculate local expansions (downwards pass)
    if (calculate_local(particles, nodes, num_particles, num_levels, num_nodes, P) != 0) {
        fprintf(stderr, "Error calculating local expansions\n");
        free(particles);
        for (int i = 0; i < num_nodes; i++) {
            free(nodes[i].expansions);
        }
        free(nodes);
        return 1;
    }

    free(particles);
    for (int i = 0; i < num_nodes; i++) {
        free(nodes[i].expansions);
    }
    free(nodes);
    return 0;
}