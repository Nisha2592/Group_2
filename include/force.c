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
    for (i = 0; i < natoms; ++i) {
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
    for (i = 0; i < natoms; ++i) {
        vx[i] += constant * fx[i];
        vy[i] += constant * fy[i];
        vz[i] += constant * fz[i];
    }
}



/* compute forces */
void force(mdsys_t *sys) 
{
    int i, j;
    double r, ffac;
    double rx1, ry1, rz1, rcsq, rsq, rxi, ryi, rzi;
    double c12, c6;
    double *fx = sys->fx;
    double *fy = sys->fy;
    double *fz = sys->fz;
    double *rx = sys->rx;
    double *ry = sys->ry;
    double *rz = sys->rz;
    double natoms = sys->natoms;
    double epsilon = sys->epsilon;
    double sigma = sys->sigma;
    double rcut = sys->rcut;

    /* zero energy and forces */
    sys->epot = 0.0;
    azzero(fx, natoms);
    azzero(fy, natoms);
    azzero(fz, natoms);

    /* precompute some constants */
    c12 = 4.0 * epsilon * pow(sigma, 12.0);
    c6 = 4.0 * epsilon * pow(sigma, 6.0);
    rcsq = rcut * rcut;
    double Box = sys->box;
    double halfbox = 0.5 * sys->box;

    for (i = 0; i < (natoms) - 1; ++i) {
        rxi = sys->rx[i];
        ryi = sys->ry[i];
        rzi = sys->rz[i];
        for (j = i + 1; j < (natoms); ++j) {

            /* get distance between particle i and j */
            rx1 = pbc(rxi - rx[j], halfbox, Box);
            ry1 = pbc(ryi - ry[j], halfbox, Box);
            rz1 = pbc(rzi - rz[j], halfbox, Box);
            rsq = rx1 * rx1 + ry1 * ry1 + rz1 * rz1;

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
                fx[i] += rx1 * ffac;
                fy[i] += ry1 * ffac;
                fz[i] += rz1 * ffac;
                fx[j] -= rx1 * ffac;
                fy[j] -= ry1 * ffac;
                fz[j] -= rz1 * ffac;
            }
        }
    }
}