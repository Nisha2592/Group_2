#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>
#include "ljmd.h"


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
    //force(sys);
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
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int i, j;
    double r, ffac;
    double rcsq, rsq;
    double c12, c6, ssigma12, ssigma6;
    double rx;
    double ry;
    double rz;
    double natoms = sys->natoms;
    double epsilon = sys->epsilon;
    double sigma = sys->sigma;
    double rcut = sys->rcut;

    /* zero energy and forces */
    double epot = 0.0;
    azzero(sys->cx, sys->natoms);
    azzero(sys->cy, sys->natoms);
    azzero(sys->cz, sys->natoms);

    //Broadcast particle positions to all ranks
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* precompute some constants */
    ssigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
    ssigma12 = ssigma6 * ssigma6;
    c12 = 4.0 * epsilon * ssigma12;
    c6 = 4.0 * epsilon * ssigma6;
    rcsq = rcut * rcut;
    double Box = sys->box;
    double halfbox = 0.5 * sys->box;

    int ii;

    for (i = 0; i < (sys->natoms) - 1; i+=size) {
        //round-robin distribution of data
        ii = i + rank;
        if (ii >= (sys->natoms-1)){
            //out of bound
            break;
        }

        for (j = ii + 1; j < (sys->natoms); ++j) {
            /* get distance between particle i and j */
            rx = pbc(sys->rx[ii] - sys->rx[j], halfbox, Box);
            ry = pbc(sys->ry[ii] - sys->ry[j], halfbox, Box);
            rz = pbc(sys->rz[ii] - sys->rz[j], halfbox, Box);
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
                sys->cx[ii] += rx * ffac;
                sys->cy[ii] += ry * ffac;
                sys->cz[ii] += rz * ffac;

                sys->cx[j] -= rx * ffac;
                sys->cy[j] -= ry * ffac;
                sys->cz[j] -= rz * ffac;
            }
        }
    }

    //Reduce forces from all ranks into global arrays
    MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    //Reduce potential energy across ranks
    MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}
