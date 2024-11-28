#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <sys/time.h>
#include "ljmd.h"

#if defined(_OPENMP)
#include <omp.h>
#endif


const double kboltz=0.0019872067;
const double mvsq2e=2390.05736153349;


void velverlet_step1(mdsys_t *sys) {
    int i;
    double natoms = sys->natoms;
    double mass = sys->mass;
    double dt = sys->dt;
    double *fx = sys->fx;
    double *fy = sys->fy;
    double *fz = sys->fz;
    double *vx = sys->vx;
    double *vy = sys->vy;
    double *vz = sys->vz;
    double *rx = sys->rx;
    double *ry = sys->ry;
    double *rz = sys->rz;
    double constant = 0.5 * dt * 1.0 / mvsq2e * 1.0 / mass;

    #if defined(_OPENMP)
    #pragma omp parallel for simd
    #endif
    for (i = 0; i < sys->natoms; ++i) {
        vx[i] += constant * fx[i];
        vy[i] += constant * fy[i];
        vz[i] += constant * fz[i];
        rx[i] += dt * sys->vx[i];
        ry[i] += dt * sys->vy[i];
        rz[i] += dt * sys->vz[i];
    }

}


void velverlet_step2(mdsys_t *sys) {
    /* compute forces and potential energy */
    int i;
    double natoms = sys->natoms;
    double mass = sys->mass;
    double dt = sys->dt;
    double *fx = sys->fx;
    double *fy = sys->fy;
    double *fz = sys->fz;
    double *vx = sys->vx;
    double *vy = sys->vy;
    double *vz = sys->vz;
    double constant = 0.5 * dt * 1.0 / mvsq2e * 1.0 / mass;

    #if defined(_OPENMP)
    #pragma omp parallel for simd
    #endif
    for (i = 0; i < sys->natoms; ++i) {
        vx[i] += constant * fx[i];
        vy[i] += constant * fy[i];
        vz[i] += constant * fz[i];
    }
}


/* compute forces */
void force(mdsys_t *sys) 
{
    double epot=0.0;
    int size = 1;
    int rank = 0;
    #ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    #endif
    double epsilon = sys->epsilon;
    double sigma = sys->sigma;
    double rcut = sys->rcut;

    #ifdef USE_MPI
    //Broadcast particle positions to all ranks
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif

    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:epot)
    #endif
    {
    /*declared inside the threaded region. 
    Each thread will have its own local copy. */    
    double *fx, *fy, *fz;
    int i,j;
    double r, ffac;
    double rcsq, rsq;
    double c12, c6, ssigma12, ssigma6;

    #if defined(_OPENMP)
    int tid = omp_get_thread_num();
    #else
    int tid = 0;
    #endif

    fx = sys->cx + (tid * sys->natoms);
    fy = sys->cy + (tid * sys->natoms);
    fz = sys->cz + (tid * sys->natoms);

    /* zero energy and forces */
    //double epot = 0.0;
    azzero(fx, sys->natoms);
    azzero(fy, sys->natoms);
    azzero(fz, sys->natoms);

    /* precompute some constants */
    ssigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
    ssigma12 = ssigma6 * ssigma6;
    c12 = 4.0 * epsilon * ssigma12;
    c6 = 4.0 * epsilon * ssigma6;
    rcsq = rcut * rcut;
    double Box = sys->box;
    double halfbox = 0.5 * sys->box;

    int ii;
    int stride = size * sys->nthreads;
    for (i = 0; i < (sys->natoms) - 1; i+=stride) {
        //round-robin distribution of iterations to processes
        ii = i + rank * sys->nthreads + tid; 
        if (ii >= (sys->natoms-1)){
            //out of bound
            break;
        }

        for (j = ii + 1; j < (sys->natoms); ++j) {
            /* get distance between particle ii and j */
            double rx = pbc(sys->rx[ii] - sys->rx[j], halfbox, Box);
            double ry = pbc(sys->ry[ii] - sys->ry[j], halfbox, Box);
            double rz = pbc(sys->rz[ii] - sys->rz[j], halfbox, Box);
            rsq = rx * rx + ry * ry + rz * rz;

            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
                double r2inv, r6inv, r12inv;
                
                // Precompute reciprocal of distance squared
                r2inv = 1.0 / rsq;
                r6inv = r2inv * r2inv * r2inv; // r^(-6)
                r12inv = r6inv * r6inv;        // r^(-12)
                    
                // Force factor
                ffac = (12.0 * c12 * r12inv - 6.0 * c6 * r6inv) * r2inv;

                // Accumulate potential energy
                epot += c12 * r12inv - c6 * r6inv;

                // Update forces
                fx[ii] += rx * ffac;
                fy[ii] += ry * ffac;
                fz[ii] += rz * ffac;

                fx[j] -= rx * ffac;
                fy[j] -= ry * ffac;
                fz[j] -= rz * ffac;
            }
        }
    }
    /*Wait for the threads to finish reduction the forces.*/
    #if defined (_OPENMP)
    #pragma omp barrier
    #endif

    i = (sys->natoms/sys->nthreads) + 1;
    int fromidx = tid * i;
    int toidx = fromidx + i;
    if (toidx > sys->natoms) toidx = sys->natoms;

    /*This is still in the paralllel region. Reduce forces
    from threads into the storage corresponding 
    to the master thread.*/
    for (int i = 1; i < sys->nthreads; ++i){
      int offs = i*sys->natoms;
      for(int k = fromidx; k < toidx; ++k){
        sys->cx[k] += sys->cx[offs + k];
        sys->cy[k] += sys->cy[offs + k];
        sys->cz[k] += sys->cz[offs + k];
      }
    }

    }//end openmp parallel region
    
    #ifdef USE_MPI
    //Reduce forces from all ranks into global arrays
    MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    //Reduce potential energy across ranks
    MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    #else
    for (int i = 0; i < sys->natoms; i++){
      sys->fx[i] = sys->cx[i];
      sys->fy[i] = sys->cy[i];
      sys->fz[i] = sys->cz[i];
    } 

    sys->epot = epot;

    #endif

}
