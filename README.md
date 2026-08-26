# fmm-openmp

A parallel implementation of the Fast Multipole Method, by [Clint Goodwin](https://github.com/cbgoodwin) and [Peter Scott](https://github.com/plscott29).

The FMM evaluates pairwise particle interactions

$$p_i = \sum_{j=1}^{N} K(x_i, x_j)\,q_j, \qquad K(x_i,x_j) = -\frac{\log|x_i - x_j|}{4\pi}$$

to relative precision $\epsilon$ in $O(N)$ operations rather than the $O(N^2)$ of direct evaluation, with a constant scaling as $\log(1/\epsilon)$. This implementation handles the 2D Laplace kernel in single precision, parallelized with OpenMP and tuned for NUMA hardware.

Full algorithm description, experiments, and analysis: **[report.pdf](report.pdf)**.

## Build

```sh
make all
```

Produces four executables:

| Executable | Purpose |
| --- | --- |
| `data_writer` | Generate test particle data |
| `fmm` | Run the FMM (sequential or parallel) |
| `direct` | Brute-force $O(N^2)$ evaluation, for validation and timing |
| `compare` | Relative error between two output files |

## Usage

```sh
./data_writer <num_particles> <file>
./fmm <input_file> <output_file> <num_levels> <p> <num_particles> <is_parallel> <num_threads> [first_touch]
./direct <input_file> <output_file> <num_particles> <is_parallel> <num_threads>
./compare <output_file1> <output_file2> <num_particles>
```

**`fmm` parameters.** `num_levels` sets tree depth; `p` is the length of each multipole and local expansion; `is_parallel` selects the parallel path and, when `1`, uses `num_threads`; `first_touch` optionally disables NUMA first-touch placement (defaults to `1`, turn off with `0`).

`compare` prints the relative error of the second file against the first, normalized by total charge across all particles.

### End-to-end

```sh
make all
./data_writer 1000 data.csv
./fmm data.csv fmm_result.csv 4 7 1000 1 4
> Total time elapsed: 0.00906
./direct data.csv direct_result.csv 1000 0 1
> Total time elapsed: 0.04168
./compare direct_result.csv fmm_result.csv 1000
> Relative precision epsilon between direct_result.csv and fmm_result.csv: 5.13874e-06
```

## Approach

The algorithm runs in four largely sequential stages, each parallel internally:

1. **Tree construction** — recursively partition the domain into quadrants down to a fixed depth `num_levels`, giving a uniform tree with $4^l$ leaves. Particles are sorted in place by an unstable four-way partition, so each node stores only an index range. Quadrants at a given depth are independent, so this parallelizes completely.
2. **Multipole expansions (upward pass)** — *outgoing from sources* across leaf nodes, then *outgoing from outgoing* level by level. Nodes within a level are independent; levels are separated by a barrier.
3. **Local expansions (downward pass)** — *incoming from outgoing* over each node's interaction list, plus *incoming from incoming* inherited from the parent. Same structure: parallel within a level, barrier between levels.
4. **Potential evaluation** — per leaf node, combining *target from incoming* with direct near-field summation over the node and its neighbors. The per-particle near-field loop stays sequential; at particle-level granularity the overhead isn't worth it.

Most of this falls out of well-placed `#pragma omp parallel for` directives, given a sequential implementation structured to expose it.

Two further optimizations:

**Sequential truncation.** Near the top of the tree there's too little work to parallelize — level 2 has only 16 nodes, small enough to sit in L1/L2. Each level requests $4^l / 128$ threads, capped at `omp_get_max_threads()`; when that falls to one, the level's node loop runs sequentially.

**First-touch placement.** Tree construction sorts the particle array in place, and nodes hold index ranges into it. Mirroring that access pattern during initialization puts each node's expansion data — and its slice of the particle array — in memory local to the thread that will compute on it.

**Data layout.** Particles are an array of structs, sorted in place. The tree is an array of node structs (level, midpoint, particle range, expansion pointers), with all expansion memory allocated as one contiguous block at startup so that neighboring nodes have neighboring expansions. A linked-list tree would have been easier to write and scaled badly.

## Results

Measured on NYU CIMS `crunchy2`/`crunchy3`: 64 cores across 4 sockets and 8 NUMA nodes, 6 MiB L3 shared per 8-core die.

**FMM vs. direct evaluation** ($P = 8$, 16 threads, $l$ chosen for near-optimal $B$):

| $N$ | $l$ | Direct serial | Direct parallel | FMM serial | FMM parallel |
| --- | --- | --- | --- | --- | --- |
| 2,048 | 3 | 2.370 s | 0.338 s | 0.296 s | 0.126 s |
| 8,192 | 4 | 197.105 s | 5.718 s | 2.510 s | 0.653 s |
| 32,768 | 5 | — | 100.114 s | 6.346 s | 1.983 s |
| 131,072 | 6 | — | 650.621 s | 15.412 s | 5.025 s |

Both FMM variants scale as $O(N)$ against the direct method's $O(N^2)$, at matching accuracy ($\epsilon \approx 10^{-7}$, single precision). Parallel speedup fell short of expectations — contention on shared machines plus the genuinely sequential portions of the algorithm.

**Speedup vs. thread count** (row 1 is baseline wall time, not speedup):

| threads | $N$ = 2,048 | $N$ = 8,192 | $N$ = 32,768 | $N$ = 131,072 |
| --- | --- | --- | --- | --- |
| 1 | (0.142 s) | (0.626 s) | (3.938 s) | (10.836 s) |
| 2 | 1.667 | 1.883 | 1.645 | 1.944 |
| 4 | 1.332 | 2.409 | 2.960 | 3.682 |
| 8 | 0.317 | 3.217 | 3.463 | 6.644 |
| 16 | 0.113 | 3.816 | 5.566 | 7.903 |

Speedup improves with problem size at every thread count, and marginal thread efficiency falls off approaching 16. At $N = 2{,}048$ the parallel version is *slower* than serial past 4 threads — thread management cost can't amortize over a sub-second run.

**Optimal leaf density** $B$ (particles per leaf box, $P$ set for $\epsilon \approx 10^{-6}$):

| $B$ | $N$ = 8,192 | $N$ = 32,768 | $N$ = 131,072 |
| --- | --- | --- | --- |
| 0.125 | 21.920 s | 39.933 s | 91.410 s |
| 0.5 | 5.381 s | 17.803 s | 7.018 s |
| 2 | 3.020 s | 4.490 s | 1.686 s |
| 8 | 0.499 s | 1.988 s | 1.360 s |
| **32** | **0.131 s** | **0.559 s** | **1.021 s** |
| 128 | 0.171 s | 2.407 s | 7.505 s |
| 512 | 1.648 s | 10.092 s | 9.544 s |
| 2,048 | 19.191 s | 57.202 s | 170.288 s |

Every problem size optimizes at $B = 32$, confirming Greengard and Gropp's $B_{\text{opt}} \approx 30$ for single precision. $B$ trades direct near-field evaluations against expansion computations: at $B = N$ the FMM degenerates to the direct $O(N^2)$ calculation, while at $B = 0.125$ most particles sit alone in a box and many boxes are empty, so expansion work is wasted on nothing.

**First-touch placement** (1,000,000 particles, $l = 8$, $P = 10$):

| threads | No first-touch | First-touch |
| --- | --- | --- |
| 4 | 31.013 s | 34.557 s |
| 8 | 16.911 s | 18.878 s |
| 16 | 10.632 s | **9.440 s** |
| 32 | 6.280 s | **5.738 s** |

First-touch only pays off at 16 threads and beyond. A NUMA node on `crunchy2` has 8 cores, so at 8 threads or fewer everything fits on one node and the extra initialization — copying the particle array, re-initializing expansion memory — is pure overhead. Past that, threads spill across nodes and locality starts to matter.

## Takeaways

- **Algorithm choice sets the ceiling on speedup.** An algorithm that parallelizes naturally beats a heavily tuned naive one with far less effort.
- **Parameters matter as much as parallelism.** A badly chosen $B$ produces an FMM that loses to a thoughtfully parallelized brute-force evaluation.
- **Locality optimizations aren't free.** First-touch costs more than it returns below the NUMA-node threshold; knowing when an optimization pays is part of the work.

## References

1. L. Greengard and V. Rokhlin, [A fast algorithm for particle simulations](https://doi.org/10.1016/0021-9991(87)90140-9), *J. Comput. Phys.* **73**(2), 1987.
2. L. Greengard and W. D. Gropp, [A parallel version of the fast multipole method](https://doi.org/10.1016/0898-1221(90)90349-O), *Computers & Mathematics with Applications* **20**(7), 1990.
3. N. Engheta, W. D. Murphy, V. Rokhlin, M. S. Vassiliou, [The fast multipole method (FMM) for electromagnetic scattering problems](https://doi.org/10.1109/8.144597), *IEEE Trans. Antennas Propag.* **40**(6), 1992.
4. P. G. Martinsson, *Fast Direct Solvers for Elliptic PDEs*, SIAM, 2020.
5. B. Zhang, J. Huang, N. P. Pitsianis, X. Sun, [RECFMM: Recursive parallelization of the adaptive fast multipole method](https://doi.org/10.4208/cicp.OA-2015-0005), *Commun. Comput. Phys.* **20**(2), 2016.

---

Clint Goodwin and Peter Scott · April 2026 · CSCI-GA 3033: Multicore Processors
