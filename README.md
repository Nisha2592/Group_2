# Lennard-Jones Molecular Dynamic Simulator
This work presents a optimized, and parallelized Lennard-Jones Molecular Dynamics simulator. The code is designed to efficiently simulate large-scale molecular systems.
## Collaborators:

        | Christian Tica         --> openMP                ---> GitHub: cristiano-mhpc

        | Nicola Aladrah         --> Optimization          ---> GitHub: NicolaAdrah

        | Nisha Nisha            --> MPI                   ---> GitHub: Nisha2592

## Compiler Options:
This code offers various compilation options, including serial, optimized, MPI, OpenMP, and a hybrid OpenMP-MPI approach. The provided `CMakeLists.txt` file facilitates the build process by using:
        * `cmake -S . -B build`  
        * `make serial`         --> To compile the serial code  
        * `make mpi`            --> To compile the MPI version  
        * `make openmp`         --> To compile the openMP version  
        * `make hybrid`         --> To compile the Hybrid MPI+openMP version  
Note that both the OpenMP and MPI versions are built upon the optimized codebase. To enable unit testing for force, utilities, and Verlet functions, add `-DENABLE_TESTING=ON` to the CMake cache.  
To optimize performance across all implementations, we utilized the following compiler flags:  

        `-Wall -O3 -ffast-math -fexpensive-optimizations -ffp-contract=fast -msse3`  

These flags enable aggressive optimization techniques, including loop unrolling, function inlining, and vectorization.

## Benchmark:
### Optimization stage:
To evaluate the performance impact of optimization, we employed Valgrind and Kcachegrind to profile the optimized and non-optimized code. The profiling process involved running the following:  
* `valgrind --tool=callgrind ./our_executable`
* `kcachegrind callgrind.out.<pid>`

In the case of non-optimized code, we have the following table and call graph,  

        | function |    calls    |
        |:--------:|:-----------:|
        |   pbc    | 346'714'668 |
        |   pow    | 237'125'560 |
        |   sqrt   | 115'571'556 |

<div style="text-align: center;">
    <img src="benchmark_opt/non_opt.jpg" alt="non-opt" width="1000"/>
</div>


By applying Newton's Third Law, the `pbc()` function becomes solely dependent on `j`, as the `i`-dependent calculations are precomputed before entering the `j` loop. Thus, studying argon for 108 atoms, we reduced the timing from 300.780 s to 146 s. We employ strength reduction to replace costly pow() function calls with more efficient multiplication operations. Additionally, we avoid unnecessary sqrt() calculations by comparing squared values directly within conditional statements. Adding the previous optimizations together reduced the timing from 300.780 s to 80 s. We can see this easily from the following call graph,

<div style="text-align: center;">
    <img src="benchmark_opt/opt.jpg" alt="non-opt" width="1000"/>
</div>

where now the calls for pbc is negligible and the force function only call itself. We plot the timing for both the optimized and non-optimized codes,

<div style="text-align: center;">
    <img src="benchmark_opt/benchmark.png" alt="benchmark" width="1000"/>
</div>

While optimizing larger systems can significantly boost their performance, smaller systems might not see much improvement.  

## MPI

## OpenMP

## Hybrid