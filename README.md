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
This report analyzes the performance of an optimized computational code executed on varying numbers of atoms using different configurations: An optimized code without message passing interface and with MPI employing 2, 4, 6, 8 and 10 processors. In this implementation, MPI is primarily leveraged to accelerate the computation of forces by distributing the workload across multiple processors. The Force function is adapted to handle this parallelization through the following steps:
* Atom positions are shared across all processes using broadcast communication.
* Process-specific indices are defined to control loop iterations, ensuring tasks are divided based on the total number of processes and their individual ranks.
* Dedicated buffers are utilized to store the computed forces for each process.
* Finally, the computed forces and total potential energy from all processes are gathered and consolidated into the rank 0 process using reduction operations.

This approach is layered on top of the pre-existing optimized codebase to achieve enhanced performance and scalability. The benchmarking data and plot is structured as follows:

| Number of Atoms |   W/O MPI   |   MPI (2 Procs)   |   MPI (4 Procs)   |   MPI (6 Procs)   |   MPI (8 Procs)   |   MPI (10 Procs)|
|:---------------:|:-----------:|:----------------:|:----------------:|:----------------:|:----------------:|:------------------:|
|       108       |    1.610    |       0.621      |       0.805      |       0.616      |      0.5118      |       0.478       |
|      2916       |   66.781    |      33.94       |      35.297      |      24.049      |      18.182      |      14.986       |
|     78732       |   600.235   |     303.002      |      342.45      |      257.68      |     211.958      |      176.681      |

<div style="text-align: center;">
      <img src="benchmark_mpi/benchmark.png" alt="benchmark" width="1000"/>
  </div>

### Analysis on basis of scalability:
* Larger systems (e.g. 78732 atoms) benefit more significantly from MPI parallelization due to the increased workload justifying the overhead of communication between processors. 
* Smaller systems (e.g. 108 atoms) exhibit diminishing returns when scaling up the number of processors.

### Analysis on basis of efficiency:
* At 108 atoms, the time reduction between 8 and 10 processors is marginal, which indicates a point of saturation where additional processors contribute little to performace.
* For 78732 atoms, the improvement remains notable even with 10 processors, which suggest that workload is large enough to benefit from additional parallel resources.

## OpenMP

## Hybrid
