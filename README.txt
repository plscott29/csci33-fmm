Parallel implementation of Fast Multiple Method (FMM)

To compile via the Makefile, run:

make all

This will generate four executables:
    1. fmm
    2. direct
    3. compare
    4. data_writer

The data_writer exectuble will generate a file with test particle data, with intended usage:

./data_writer <num_particles> <file>

The fmm exectuble runs the FMM algorithm on a set of input particles, with intended usage:

./fmm <input_file> <output_file> <num_levels> <p> <num_particles> <is_parallel> <num_threads> [first_touch]

The input_file parameter is a filepath to the input particle data file (e.g., generated via data_writer).
The output_file parameter is the filepath to where the output particle data will be written after completion.
The num_levels parameter specifies the level of depth of the tree.
The p parameter specifies the number of multipole and local expansions (each of length p) to calculate in the algorithm.
The num_particles parameter specifies how many particles in input_file to evaluate.
The is_parallel parameter specifies whether to run the parallel implementation of the FMM algorithm,
and if it is 1 it will use the num_threads parameter to specify how many threads to run.
The first_touch parameter is an optional parameter that specifies whether to use NUMA first-touch data placement
in the parallel FMM implementation (defaults to 1; turn off with 0).

The direct executable runs a direct pair-wise potential evaluation algorithm on a set of input particles, with intended usage:

./direct <input_file> <output_file> <num_particles> <is_parallel> <num_threads>

The input_file parameter is a filepath to the input particle data file (e.g., generated via data_writer).
The output_file parameter is the filepath to where the output particle data will be written after completion.
The num_particles parameter specifies how many particles in input_file to evaluate.
The is_parallel parameter specifies whether to run the parallel implementation of the FMM algorithm,
and if it is 1 it will use the num_threads parameter to specify how many threads to run.

The compare executable can be used to compare the results of two separate runs of potential evaluation algorithms,
(e.g., to validate the results of FMM against the brute-force pair-wise potential calculation) with intended usage:

./compare <output_file1> <output_file2> <num_particles>

The output_file1 parameter is the filepath to the first output particle data to read.
The output_file2 parameter is the filepath to the second output particle data to read.
The num_particles parameter specifies how many particles in input_file to evaluate.

The compare executable will print out the relative error of the particles in output_file2
compared to the partiles in output_file1, normalized by the total charge across all particles.

An end-to-end run would take the form of (1) generate test data, (2) run FMM to evaluate
particle potentials, (3) run direct evaluation of particle potentials, and (4) compare results.
This would look like:

make all
./data_writer 1000 data.csv
./fmm data.csv fmm_result.csv 4 7 1000 1 4
> Total time elapsed: 0.00906
./direct data.csv direct_result.csv 1000 0 1
> Total time elapsed: 0.04168
./compare direct_result.csv fmm_result.csv 1000
Relative precision epslon between direct_result.csv and fmm_result.csv: 5.13874e-06