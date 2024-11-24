#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
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



/* compute forces */
// void force(mdsys_t *sys)
// {
//     double r,ffac, c12, c6;
//     double rx,ry,rz, rcsq, rsq;
//     int i,j;

//     /* zero energy and forces */
//     sys->epot=0.0;
//     azzero(sys->fx,sys->natoms);
//     azzero(sys->fy,sys->natoms);
//     azzero(sys->fz,sys->natoms);

//     // --------------- OPTIMIZATION ATTEMPT II ---------------
//     // Removing expensive math using
//     c12 = 4.0 * sys->epsilon * pow(sys->sigma,12.0);
//     c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
//     rcsq = sys->rcut * sys->rcut;

//     // --------------- OPTIMIZATION ATTEMPT I -----------------
//     // Newton 3rd law
//     for(i=0; i < (sys->natoms) - 1; ++i) {
//         for(j= i + 1; j < (sys->natoms); ++j) {

//             /* particles have no interactions with themselves */
//             // if (i==j) continue;

//             /* get distance between particle i and j */
//             rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
//             ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
//             rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
//             rsq = rx*rx + ry*ry + rz*rz;


//             /* compute force and energy if within cutoff */
//             if (rsq < rcsq) {
//                 double r6, rinv;
//                 rinv = 1.0/rsq;
//                 r6 = rinv * rinv * rinv;
//                 ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;
//                 sys->epot += r6 * (c12 * r6 - c6);

//                 // Consider action in i direction,
//                 sys->fx[i] += rx*ffac;
//                 sys->fy[i] += ry*ffac;
//                 sys->fz[i] += rz*ffac;
//                 // then the reaction in the j direction.
//                 sys->fx[i] -= rx*ffac;
//                 sys->fy[i] -= ry*ffac;
//                 sys->fz[i] -= rz*ffac;
//             }
//         }
//     }
// }


void force(mdsys_t *sys)
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < (sys->natoms); ++i) {
        for(j=0; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);

                sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));

                sys->fx[i] += rx/r*ffac;
                sys->fy[i] += ry/r*ffac;
                sys->fz[i] += rz/r*ffac;
            }
        }
    }
}