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

typedef enum {
    NORTH,
    SOUTH,
    EAST,
    WEST
} DIRECTION;


/* Helper functions:
    - child_idx: given a parent node index, returns the index of the first child in the quadtree
    - parent_idx: given a child node index, returns the index of the parent node in the quadtree
    - level_start_idx: given a quadtree level, returns the index of the first node at that level
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
    while (low <= high) {
        if (particles[low].x <= node->x_mid) //if (creal(particles[low].z) <= creal(node->z_mid)) TODO: Remove x, y fields (only use z)
            low++;
        else if (particles[high].x > node->x_mid) //else if (creal(particles[high].z) > creal(node->z_mid)) TODO: Remove x, y fields (only use z)
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
    while (low <= high) {
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
    while (low <= high) {
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

// returns the index of the neighbor of the node in the provided direction, -1 if it does not exist
int search_for_neighbor(Node *nodes, int node_idx, DIRECTION direction, float distance) {
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
        case NORTH:
            while (test_idx <= idx_upper_bound) {
                jump_idx = pow(4, iteration)-level_start_idx(iteration); // +4-1=3, +16-4-1=11, +64-16-4-1=43, etc.
                test_idx = node_idx + jump_idx; // we want to see if this node is the north neighbor
                test_y_dist = cimag(nodes[test_idx].z_mid - nodes[node_idx].z_mid);
                ratio = test_y_dist/distance;
                aligned = creal(nodes[test_idx].z_mid) == creal(nodes[node_idx].z_mid);
                if (0.8 < ratio && ratio < 1.2 && aligned) {return test_idx;}
                iteration++;
            }
            break;

        case SOUTH:
            while (test_idx >= idx_lower_bound) {
                jump_idx = pow(4, iteration)-level_start_idx(iteration);
                test_idx = node_idx - jump_idx;
                test_y_dist = cimag(nodes[node_idx].z_mid - nodes[test_idx].z_mid);
                ratio = test_y_dist/distance;
                aligned = creal(nodes[test_idx].z_mid) == creal(nodes[node_idx].z_mid);
                if (0.95 < ratio && ratio < 1.05 && aligned) {return test_idx;}
                iteration++;
            }
            break;

        case EAST:
            while (test_idx <= idx_upper_bound) {
                jump_idx = 2*(pow(4, iteration)-level_start_idx(iteration));
                test_idx = node_idx + jump_idx;
                test_x_dist = creal(nodes[test_idx].z_mid - nodes[node_idx].z_mid);
                ratio = test_x_dist/distance;
                aligned = cimag(nodes[test_idx].z_mid) == cimag(nodes[node_idx].z_mid);
                if (0.8 < ratio && ratio < 1.2 && aligned) {return test_idx;}
                iteration++;
            }
            break;

        case WEST:
            while (test_idx >= idx_lower_bound) {
                jump_idx = 2*(pow(4, iteration)-level_start_idx(iteration));
                test_idx = node_idx - jump_idx;
                test_x_dist = creal(nodes[node_idx].z_mid - nodes[test_idx].z_mid);
                ratio = test_x_dist/distance;
                aligned = cimag(nodes[test_idx].z_mid) == cimag(nodes[node_idx].z_mid);
                if (0.95 < ratio && ratio < 1.05 && aligned) {return test_idx;}
                iteration++;
            }
            break;
    }
    return -1;
}

// returns Neighborhood struct with indices of neighboring nodes, -1 index for the node itself and non-existent neighbors
Neighborhood get_neighborhood(Node *nodes, int node_idx) {
    int sw_child_idx = child_idx(parent_idx(node_idx));
    float mid2mid = creal(nodes[sw_child_idx+2].z_mid) - creal(nodes[sw_child_idx].z_mid);
    int north, south, east, west, ne, se, nw, sw;

    switch (node_idx-sw_child_idx) {
        case 0: //search for south, west, nw, sw, se
            north = node_idx+1; east = node_idx+2; ne = node_idx+3;
            south = search_for_neighbor(nodes, node_idx, SOUTH, mid2mid);
            west = search_for_neighbor(nodes, node_idx, WEST, mid2mid);

            if (south == -1 && west == -1) {nw = -1; sw = -1; se = -1;}
            else if (south == -1) {nw = west+1; sw = -1; se = -1;}
            else if (west == -1) {nw = -1; sw = -1; se = south+2;}
            else {nw = west+1; sw = search_for_neighbor(nodes, west, SOUTH, mid2mid); se = south+2;}
            break;

        case 1: //search for north, west, ne, nw, sw
            south = node_idx-1; east = node_idx+2; se = node_idx+1;
            north = search_for_neighbor(nodes, node_idx, NORTH, mid2mid);
            west = search_for_neighbor(nodes, node_idx, WEST, mid2mid);

            if (north == -1 && west == -1) {ne = -1; nw = -1; sw = -1;}
            else if (north == -1) {ne = -1; nw = -1; sw = west-1;}
            else if (west == -1) {ne = north+2; nw = -1; sw = -1;}
            else {ne = north+2; nw = search_for_neighbor(nodes, west, NORTH, mid2mid); sw = west-1;}
            break;

        case 2: //search south, east, sw, se, ne
            north = node_idx+1; west = node_idx-2; nw = node_idx - 1;
            south = search_for_neighbor(nodes, node_idx, SOUTH, mid2mid);
            east = search_for_neighbor(nodes, node_idx, EAST, mid2mid);

            if (south == -1 && east == -1) {sw = -1; se = -1; ne = -1;}
            else if (south == -1) {sw = -1; se = -1; ne = east+1;}
            else if (east == -1) {sw = south-2; se = -1; ne = -1;}
            else {sw = south-2; se = search_for_neighbor(nodes, east, SOUTH, mid2mid); ne = east+1;}
            break;

        case 3: //search north, east, se, ne, nw
            south = node_idx-1; west = node_idx-2; sw = node_idx-3;
            north = search_for_neighbor(nodes, node_idx, NORTH, mid2mid);
            east = search_for_neighbor(nodes, node_idx, EAST, mid2mid);

            if (north == -1 && east == -1) {se = -1; ne = -1; nw = -1;}
            else if (north == -1) {se = east-1; ne = -1; nw = -1;}
            else if (east == -1) {se = -1; ne = -1; nw = north-2;}
            else {se = east-1; ne = search_for_neighbor(nodes, east, NORTH, mid2mid); nw = north-2;}
            break;

        default:
            north = -1; south = -1; east = -1; west = -1;
            ne = -1; se = -1; nw = -1; sw = -1;
    }
    return (Neighborhood) {
        .nodes = {
            {nw, north, ne},
            {west, -1, east},
            {sw, south, se}}
    };
}

// returns up to 27 well separated node indices, padded with -1 if there are less than 27 well separated indices
WellSeparated get_well_separated(Node *nodes, int node_idx) {
    Neighborhood parent_neighborhood = get_neighborhood(nodes, parent_idx(node_idx));
    WellSeparated *well_separated = malloc(sizeof(WellSeparated));
    int ws_idx = 0;

    for (int i = 0; i < 3; i++) { for (int j = 0; j < 3; j++) {
        int c_idx;
        int p_idx = parent_neighborhood.nodes[i][j];
        if (p_idx != -1) { // case when
            c_idx = child_idx(p_idx);
            for (int k = 0; k < 4; k++) {
                float distance = cabsf(nodes[node_idx].z_mid - nodes[c_idx].z_mid);
                float distance_ratio = distance*pow(2.0, nodes[node_idx].level);
                if (distance_ratio > 1.9) {well_separated->nodes[ws_idx] = c_idx; ws_idx++;}
                c_idx++;
            }
        }
    }}

    for (int i = ws_idx; i < 27; i++) {well_separated->nodes[i] = -1;}
    return *well_separated;
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

    // TODO: use first touch placement to copy particles into new array such that each leaf node's particles
    // are stored in memory close to where threads would access them in later computations (NUMA)

    // TODO: use first touch placement to initialize the expansions array chunks for each node such that
    // they are stored in memory close to where threads would access them in later computations (NUMA)

    return 0;
}


int construct_tree_parallel(Particle *particles, Node *nodes, int num_particles, int num_levels)
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
        int level_start = level_start_idx(level);
        # pragma omp parallel for private(c_idx, p_idx)
        for (int i = 0; i < num_boxes; i++) {
            p_idx = level_start + i;
            // sort particles and get the partitions of subarray of particles for this parent node
            Partition partitions = four_sort(particles, &nodes[p_idx]);

            // set four nodes corresponding to the four quadrants of this parent node
            c_idx = child_idx(p_idx);
            //printf("Level: %d, Parent: %d, Child: %d, Thread: %d\n", level, p_idx, c_idx, omp_get_thread_num());
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
        node->expansions[0] += (float) particles[j].q;

        float complex offset = particles[j].z - node->z_mid;
        for (int p = 1; p < P; p++) {
            node->expansions[p] -= cpow(offset, p) * (float) particles[j].q / (float) p;
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

        for (int r = 0; r < P; r++) {
            // accumulate powers of offset for use in binomial expansion
            float complex offset_power_acc = 1;
            
            // note: need to convert between form of outgoing expansionsin equation 7.2 and form
            // expected in Theorem 7.1 where the outgoing expansions in equation 7.2 have a 
            // leading -1/p term not included in the form expected by Theorem 7.1
            float complex q_hat_r = 0;

            // iterate backwards to fill in multipole expansions from low degree to high degree
            for (int s = r; s >= 0; s--) {
                float complex q_hat_s;
                if (s == 0)
                    q_hat_s = nodes[c_idx+i].expansions[0];
                else
                    q_hat_s = - (float) s * nodes[c_idx+i].expansions[s];

                // accumulate contribution of child multipole expansion to parent    
                q_hat_r += choose(r, s) * offset_power_acc * q_hat_s;
                offset_power_acc *= offset;
            }

            // convert back to form with -1/p factor in outgoing expansions
            // (because this form is expected by IFO per equation 6.26)
            if (r == 0)
                nodes[node_idx].expansions[0] += q_hat_r;
            else
                nodes[node_idx].expansions[r] += -q_hat_r / (float) r;
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
    for (int r = P-1; r >= 0; r--) {
        // accumulate powers of offset for use in binomial expansion
        float complex offset_power_acc = 1;

        // iterate backwards to fill in multipole expansions from low degree to high degree
        for (int p = r; p <= P-1; p++) {
            // accumulate contribution of child multipole expansion to parent
            int p_choose_r = choose(p, r); // p choose r -- TODO: use binom[(p-1)*P + (r-1)];
            incoming_expansions[r] += p_choose_r * offset_power_acc * nodes[p_idx].expansions[P+p];
            offset_power_acc *= offset;
        }
    }
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

    WellSeparated well_separated_nodes = get_well_separated(nodes, node_idx);
    int interaction_list[27];
    int interaction_list_size = 0;
    for (int i = 0; i < 27; i++) {
        if (well_separated_nodes.nodes[i] != -1) {
            interaction_list[interaction_list_size++] = well_separated_nodes.nodes[i];
        }
    }

    for (int i = 0; i < interaction_list_size; i++) {
        int c_idx = interaction_list[i];
        float complex offset = nodes[node_idx].z_mid - nodes[c_idx].z_mid;          // current - interaction list node
        float complex reverse_offset = nodes[c_idx].z_mid - nodes[node_idx].z_mid;  // interaction list node - current

        // r = 0 case
        incoming_expansions[0] += nodes[c_idx].expansions[0] * clogf(offset);
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

/*
    input:
        - particles: array of Particle structs representing the source/target particles
        - nodes: tree of Nodes representing FMM hiearchical decomposition
        - node_idx: index of node in nodes array for which target from incoming expansions
        is being computed
        - P: number of terms in multipole expansion
    output:
        - particles[p_idx].p: the potential at each particle contained by this node is updated
        to include contribution from the incoming expansions at this node
    
    Compute "target from incoming" potential for each particle at this node according to:
        v^\tau = T^tfi_\tau v^\tau
    where T^tfi_\tau is the translation operator that maps incoming expansions at node
    \tau to the target potential:
        T^tfi_\tau(i, p) = (x_i - c_\tau)^(p-1)
    and v^\tau is the incoming expansion at node \tau
*/
void tfi(Particle *particles, Node *nodes, int node_idx, int P) {
    float complex *incoming_expansions = &nodes[node_idx].expansions[P];

    for (int p_idx = nodes[node_idx].start; p_idx <= nodes[node_idx].end; p_idx++) {
        float complex offset = particles[p_idx].z - nodes[node_idx].z_mid;
        float complex offset_power_acc = 1.0 + 0.0*I;
        for (int p = 0; p < P; p++) {
            particles[p_idx].p += crealf(incoming_expansions[p] * offset_power_acc);
            offset_power_acc *= offset;
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

int calculate_multipole_parallel(Particle *particles, Node *nodes, int num_particles, int num_levels, int num_nodes, int P) {
    /*
     * Upwards pass: compute outgoing expansion of each box from leaf nodes up to the root.
     * For a leaf node, compute expansion directly particles in the node ("outgoing from sources");
     * for a parent node, use outgoing expansions of its children to ("outgoing from outgoing")
     */

    // compute multipole expansions for all leaf nodes: "outgoing from sources"
    int num_leaves = pow(4, num_levels);

    #pragma omp parallel for
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
    //for (int i = (num_nodes - num_leaves) - 1; i > 0; i--) {
    //    ofo(nodes, i, P, binom);
    //}

    for (int l = num_levels-1; l > 0; l--) {
        int start_idx = level_start_idx(l); int stop_idx = level_start_idx(l+1);
        if (l >= 5) {
            #pragma omp parallel for
            for (int i = start_idx; i < stop_idx; i++) {
                ofo(nodes, i, P, binom);
            }
        }
        else {
            for (int i = start_idx; i < stop_idx; i++) {
                ofo(nodes, i, P, binom);
            }
        }
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
        memset(&nodes[i].expansions[P], 0, P*sizeof(float complex));
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

int calculate_local_parallel(Particle *particles, Node *nodes, int num_particles, int num_levels, int num_nodes, int P) {
    /*
     * Downwards pass: compute incoming expansion for every node in a pass over all nodes,
     * going from larger->smaller. For each node, combine the incoming expansion of its parent
     * with the contributions from the outgoing expansions of all boxes in its interaction list.
    */

    // set incoming expansion of root & level 1 nodes to zero
    for (int i = 0; i < (1 + 4); i++) { // root node + 4 level 1 nodes
        memset(&nodes[i].expansions[P], 0, P*sizeof(float complex));
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

        // TODO: only run this in parallel above a certain threshold of nodes in level
        // since there is overhead to parallelization and levels with very few nodes
        // may not benefit from parallelization

        #pragma omp parallel for
        for (int i = 0; i < nodes_in_level; i++) {
            int node_idx = level_start_idx + i;
            ifi(nodes, node_idx, P, binom);
            ifo(nodes, node_idx, P);
        }
        // implicit barrier: all threads must finish processing level l before level l+1
    }

    free(binom);
    return 0;
}

int evaluate_potentials(Particle *particles, Node *nodes, int num_particles, int num_levels, int num_nodes, int P) {
    /*
     * Evaluate potentials at each particle by combining contributions from
     * local expansions at the leaf node containing the particle and
     * direct interactions with nearby particles in the node itself and neighboring nodes
     * not captured by the multipole expansions.
     */

    int num_leaves = pow(4, num_levels);
    for (int i = 0; i < num_leaves; i++) {
        int leaf_idx = num_nodes - num_leaves + i;
        Node *node = &nodes[leaf_idx];
        Neighborhood neighborhood = get_neighborhood(nodes, leaf_idx);

        // initialize potentials at each particle in this leaf node to zero
        for (int p_idx = node->start; p_idx <= node->end; p_idx++) {
            particles[p_idx].p = 0.0f;
        }

        // compute targets from incoming expansions for each particle in this leaf node
        tfi(particles, nodes, leaf_idx, P);

        // compute near-field contributions from particles in this leaf node and neighboring
        // nodes via direct contribution
        for (int p_idx = node->start; p_idx <= node->end; p_idx++) {
            Particle *p = &particles[p_idx];

            // evaluate contributions from neighboring nodes
            for (int x = 0; x < 3; x++) {
                for (int y = 0; y < 3; y++) {
                    if (neighborhood.nodes[x][y] != -1) {
                        int neighbor_idx = neighborhood.nodes[x][y];
                        Node *neighbor = &nodes[neighbor_idx];

                        for (int q_idx = neighbor->start; q_idx <= neighbor->end; q_idx++) {
                            Particle q = particles[q_idx];
                            float complex offset = p->z - q.z;
                            p->p += crealf(q.q * clogf(offset));
                        }
                    }
                }
            }

            // evaluate contributions from particles in this node
            for (int q_idx = node->start; q_idx <= node->end; q_idx++) {
                Particle q = particles[q_idx];
                if (q_idx != p_idx) {
                    float complex offset = p->z - q.z;
                    p->p += crealf(q.q * clogf(offset));
                }
            }
        }
    }

    return 0;
}

int evaluate_potentials_parallel(Particle *particles, Node *nodes, int num_particles, int num_levels, int num_nodes, int P) {
    /*
     * Evaluate potentials at each particle by combining contributions from
     * local expansions at the leaf node containing the particle and
     * direct interactions with nearby particles in the node itself and neighboring nodes
     * not captured by the multipole expansions.
     * 
     * The parallel implementation should be parallelized over leaf nodes, which should be perfectly parallel
     * since there are no dependencies between leaf nods at this stage of the algorithm.
     */

    int num_leaves = pow(4, num_levels);
    #pragma omp parallel for
    for (int i = 0; i < num_leaves; i++) {
        int leaf_idx = num_nodes - num_leaves + i;
        Node *node = &nodes[leaf_idx];
        Neighborhood neighborhood = get_neighborhood(nodes, leaf_idx);

        // initialize potentials at each particle in this leaf node to zero
        for (int p_idx = node->start; p_idx <= node->end; p_idx++) {
            particles[p_idx].p = 0.0f;
        }

        // step 5: far-field potential via "targets from incoming" expansions for each particle
        tfi(particles, nodes, leaf_idx, P);

        // step 6: near-field contributions from particles in neighbors (including this node) via direct contribution
        for (int p_idx = node->start; p_idx <= node->end; p_idx++) {
            Particle *p = &particles[p_idx];

            // evaluate contributions from neighboring nodes
            for (int x = 0; x < 3; x++) {
                for (int y = 0; y < 3; y++) {
                    if (neighborhood.nodes[x][y] != -1) {
                        int neighbor_idx = neighborhood.nodes[x][y];
                        Node *neighbor = &nodes[neighbor_idx];

                        for (int q_idx = neighbor->start; q_idx <= neighbor->end; q_idx++) {
                            Particle q = particles[q_idx];
                            float complex offset = p->z - q.z;
                            p->p += crealf(q.q * clogf(offset));
                        }
                    }
                }
            }

            // evaluate contributions from particles in this node
            for (int q_idx = node->start; q_idx <= node->end; q_idx++) {
                Particle q = particles[q_idx];
                if (q_idx != p_idx) {
                    float complex offset = p->z - q.z;
                    p->p += crealf(q.q * clogf(offset));
                }
            }
        }
    }

    return 0;
}

int brute_force(Particle *particles, int num_particles) {
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
    ./fmm <input_file> <output_file> <num_levels> <p> <num_particles> <is_parallel> <num_threads>
*/
int main(int argc, char * argv[])
{
    double t_0 = omp_get_wtime();
    double t_prev;
    double t_next;
    // read inputs from command line
    if (argc != 8) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> <num_levels> <p> <num_particles> <is_parallel> <num_threads>\n", argv[0]);
        fprintf(stderr, "input_file: path to the input file containing particle data\n");
        fprintf(stderr, "output_file: path to the output file where results will be written\n");
        fprintf(stderr, "num_levels: depth of FMM expansion tree\n");
        fprintf(stderr, "P: order of multipole expansion\n");
        fprintf(stderr, "num_particles: number of particles\n");
        fprintf(stderr, "is_parallel: whether to run in parallel (1) or sequential (0)\n");
        fprintf(stderr, "num_threads: number of threads to use\n");
        return 1;
    }

    char *input_file = argv[1];
    char *output_file = argv[2];
    int num_levels = atoi(argv[3]);
    int P = atoi(argv[4]);
    int num_particles = atoi(argv[5]);
    int is_parallel = atoi(argv[6]);
    int num_threads = atoi(argv[7]);

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
    // TODO: consolidate allocations for expansions into a single large block to improve memory locality
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

    if (is_parallel) { // parallel execution
        // set # of threads for parallel execution
        omp_set_num_threads(num_threads);
        printf("Running parallel version\n");
        t_prev = omp_get_wtime();
        printf("Starting tree construction\n");
        // run step 1: tree construction & sorting
        if (construct_tree_parallel(particles, nodes, num_particles, num_levels) != 0) {
        //if (construct_tree(particles, nodes, num_particles, num_levels) != 0) {
            fprintf(stderr, "Error constructing tree\n");
            free(particles);
            for (int i = 0; i < num_nodes; i++) {
                free(nodes[i].expansions);
            }
            free(nodes);
            return 1;
        }
        t_next = omp_get_wtime();
        printf("Finished, time elapsed: %.5f\n", t_next-t_prev);

        t_prev = omp_get_wtime();
        printf("Starting multipole\n");
        // run step 2: calculate multipole expansions (upwards pass)
        if (calculate_multipole_parallel(particles, nodes, num_particles, num_levels, num_nodes, P) != 0) {
            fprintf(stderr, "Error calculating multipole expansions\n");
            free(particles);
            for (int i = 0; i < num_nodes; i++) {
                free(nodes[i].expansions);
            }
            free(nodes);
            return 1;
        }
        t_next = omp_get_wtime();
        printf("Finished, time elapsed: %.5f\n", t_next-t_prev);

        t_prev = omp_get_wtime();
        printf("Starting local\n");
        // run step 3: calculate local expansions (downwards pass)
        if (calculate_local_parallel(particles, nodes, num_particles, num_levels, num_nodes, P) != 0) {
            fprintf(stderr, "Error calculating local expansions\n");
            free(particles);
            for (int i = 0; i < num_nodes; i++) {
                free(nodes[i].expansions);
            }
            free(nodes);
            return 1;
        }
        t_next = omp_get_wtime();
        printf("Finished, time elapsed: %.5f\n", t_next-t_prev);


        t_prev = omp_get_wtime();
        printf("Starting particle potentials\n");
        // run step 4: evaluate potentials at every leaf node
        if (evaluate_potentials_parallel(particles, nodes, num_particles, num_levels, num_nodes, P) != 0) {
            fprintf(stderr, "Error evaluating potentials at leaf nodes\n");
            free(particles);
            for (int i = 0; i < num_nodes; i++) {
                free(nodes[i].expansions);
            }
            free(nodes);
            return 1;
        }
        t_next = omp_get_wtime();
        printf("Finished, time elapsed: %.5f\n", t_next-t_prev);
    } else { // sequential execution
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

        // run step 4: evaluate potentials at every leaf node
        if (evaluate_potentials(particles, nodes, num_particles, num_levels, num_nodes, P) != 0) {
            fprintf(stderr, "Error evaluating potentials at leaf nodes\n");
            free(particles);
            for (int i = 0; i < num_nodes; i++) {
                free(nodes[i].expansions);
            }
            free(nodes);
            return 1;
        }
    }
    
    // brute force calculation of particle-wise potentials (on copied particles
    // array so as not to overwrite FMM results)
    Particle * particles_bf = (Particle *)malloc(num_particles * sizeof(Particle));
    memcpy(particles_bf, particles, num_particles * sizeof(Particle));
    brute_force(particles_bf, num_particles);

    // compare brute-force and FMM results!
    float max_error = 0.0f;
    for (int i = 0; i < num_particles; i++) {
        float error = fabsf(particles[i].p - particles_bf[i].p);
        if (error > max_error) {
            max_error = error;
        }
        //printf("Particle %d:\n \tFMM potential = %f, brute-force potential = %f, error = %e\n\n",
        //       i, particles[i].p, particles_bf[i].p, fabsf(particles[i].p - particles_bf[i].p));
    }
    printf("Maximum error between FMM and brute-force potentials: %e\n", max_error);

    free(particles);
    free(particles_bf);
    for (int i = 0; i < num_nodes; i++) {
        free(nodes[i].expansions);
    }
    free(nodes);
    return 0;
}