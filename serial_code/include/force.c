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
    for (i = 0; i < sys->natoms; ++i) {
        vx[i] += constant * fx[i];
        vy[i] += constant * fy[i];
        vz[i] += constant * fz[i];
    }
}


/* compute forces */
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