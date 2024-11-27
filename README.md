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