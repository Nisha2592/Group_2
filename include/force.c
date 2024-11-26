#include <math.h>
#include <mpi.h>
#include <sys/time.h>
#include "ljmd.h"

const double kboltz=0.0019872067;
const double mvsq2e=2390.05736153349;

void velverlet_step1(mdsys_t *sys) {
    int i;
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt * sys->vx[i];
        sys->ry[i] += sys->dt * sys->vy[i];
        sys->rz[i] += sys->dt * sys->vz[i];
    }

}

void velverlet_step2(mdsys_t *sys) {
    /* compute forces and potential energy */
    //force(sys);
    int i;
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}


void force(mdsys_t *sys, int rank, int size, MPI_Comm comm) {

    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    /*set auxiliary data to 0*/
    double epot = 0.0;
    azzero(sys->cx,sys->natoms);
    azzero(sys->cy,sys->natoms);
    azzero(sys->cz,sys->natoms);

    /*broadcats particles positions from process of rank 0 to other ranks*/
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, comm);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, comm);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, comm);

    int ii;

    for(i=0; i < (sys->natoms); i+=size) {
        /*round robin distributions*/
        ii = i + rank;
        if (ii > sys->natoms) {
            break;
        }
        for(j=0; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            if (ii==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[ii] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[ii] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[ii] - sys->rz[j], 0.5*sys->box);
            // this can be avoid since we are using even powers of r
            r = sqrt(rx*rx + ry*ry + rz*rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {

                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);

                //store results to auxiliary
                epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));

                sys->cx[ii] += rx/r*ffac;
                sys->cy[ii] += ry/r*ffac;
                sys->cz[ii] += rz/r*ffac;
            }
        }
    }

    //sum forces contribution from all process to process of rank 0
    MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, comm);

    //sum potential energy contribution from all ranks to process of rank 0
    MPI_Reduce(&epot, &(sys->epot), 1, MPI_DOUBLE, MPI_SUM, 0, comm);

}
