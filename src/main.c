/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include "../include/ljmd.h"

#define BLEN 200

/* main */
int main(int argc, char **argv)
{
    int my_rank, comm_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    char ergfile[BLEN], line[BLEN];

    int nprint, i;
    
    char restfile[BLEN], trajfile[BLEN];


    FILE *fp,*traj,*erg;
    mdsys_t sys;
    double t_start;

    printf("LJMD version %3.1f\n", LJMD_VERSION);

    t_start = wallclock();

    /* read input file */
    if(!my_rank){
      if(get_a_line(stdin,line)) return 1;
      sys.natoms=atoi(line);
      if(get_a_line(stdin,line)) return 1;
      sys.mass=atof(line);
      if(get_a_line(stdin,line)) return 1;
      sys.epsilon=atof(line);
      if(get_a_line(stdin,line)) return 1;
      sys.sigma=atof(line);
      if(get_a_line(stdin,line)) return 1;
      sys.rcut=atof(line);
      if(get_a_line(stdin,line)) return 1;
      sys.box=atof(line);
      if(get_a_line(stdin,restfile)) return 1;
      if(get_a_line(stdin,trajfile)) return 1;
      if(get_a_line(stdin,ergfile)) return 1;
      if(get_a_line(stdin,line)) return 1;
      sys.nsteps=atoi(line);
      if(get_a_line(stdin,line)) return 1;
      sys.dt=atof(line);
      if(get_a_line(stdin,line)) return 1;
      nprint=atoi(line);
    }    

    //broadcast input data to other processes
    MPI_Bcast(&sys.epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.sigma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.rcut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.box, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nprint, 1, MPI_INT, 0, MPI_COMM_WORLD);


    /* allocate memory */
    allocate_mem(&sys);

    /* read restart. Only process zero.*/
    if (!my_rank){
      fp=fopen(restfile,"r");
      if(fp) {
        for (i=0; i<sys.natoms; ++i) {
          fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        for (i=0; i<sys.natoms; ++i) {
          fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(fp);
        azzero(sys.fx, sys.natoms);
        azzero(sys.fy, sys.natoms);
        azzero(sys.fz, sys.natoms);

      } else {
        perror("cannot read restart file");
        return 3;
      }

    }

    /* initialize forces and energies across all processes.*/
    sys.nfi=0;
    force(&sys, my_rank, comm_size, MPI_COMM_WORLD);
    ekin(&sys);

    if (!my_rank){
      erg=fopen(ergfile,"w");
      traj=fopen(trajfile,"w");

      printf("Startup time: %10.3fs\n", wallclock()-t_start);
      printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
      printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
      output(&sys, erg, traj);
    }

    /* reset timer */
    t_start = wallclock();

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if (!my_rank){
          if ((sys.nfi % nprint) == 0){
              output(&sys, erg, traj);
          }
        }
        /* propagate system and recompute energies */
        velverlet_step1(&sys);
        force(&sys, my_rank, comm_size, MPI_COMM_WORLD);
        velverlet_step2(&sys);
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    if(!my_rank){
      printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);
    }
    

    if(!my_rank){
      fclose(erg);
      fclose(traj);
    }

    cleanup(&sys);

    MPI_Finalize();
 
    return 0;
}
