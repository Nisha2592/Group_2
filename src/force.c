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
void force(mdsys_t *sys) 
{
    double r, ffac;
    double rx, ry, rz, rcsq, rsq;
    double c12, c6;
    int i, j;

    /* zero energy and forces */
    sys->epot = 0.0;
    azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);

    /* precompute some constants */
    c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
    c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
    rcsq = sys->rcut * sys->rcut;

    for (i = 0; i < (sys->natoms) - 1; ++i) {
        for (j = i + 1; j < (sys->natoms); ++j) {

            /* get distance between particle i and j */
            rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
            ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
            rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
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
                sys->epot += c12 * r12inv - c6 * r6inv;

                // Update forces
                sys->fx[i] += rx * ffac;
                sys->fy[i] += ry * ffac;
                sys->fz[i] += rz * ffac;
                sys->fx[j] -= rx * ffac;
                sys->fy[j] -= ry * ffac;
                sys->fz[j] -= rz * ffac;
            }
        }
    }
}