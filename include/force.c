#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include "ljmd.h"


const double kboltz=0.0019872067;
const double mvsq2e=2390.05736153349;


void velverlet_step1(mdsys_t *sys) {
    
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
    for (int i = 0; i < sys->natoms; ++i) {
        vx[i] += constant * fx[i];
        vy[i] += constant * fy[i];
        vz[i] += constant * fz[i];
        rx[i] += dt * sys->vx[i];
        ry[i] += dt * sys->vy[i];
        rz[i] += dt * sys->vz[i];
    }

}


void velverlet_step2(mdsys_t *sys) {    
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
    for (int i = 0; i < sys->natoms; ++i) {
        vx[i] += constant * fx[i];
        vy[i] += constant * fy[i];
        vz[i] += constant * fz[i];
    }
}


/* compute forces */
void force(mdsys_t *sys) 
{
    double epot = 0.0;

    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:epot)
    #endif
    {
    double *fx, *fy, *fz;
    #if defined(_OPENMP)
    int tid = omp_get_thread_num();
    #else
    int tid = 0;
    #endif

    double r, ffac;
    double rx1, ry1, rz1, rcsq, rsq;
    double c12, c6, ssigma12, ssigma6;
    double epsilon = sys->epsilon;
    double sigma = sys->sigma;
    double rcut = sys->rcut;

    fx = sys->cx + (tid * sys->natoms);
    fy = sys->cy + (tid * sys->natoms);
    fz = sys->cz + (tid * sys->natoms);

    /* zero energy and forces */
    azzero(fx, sys->natoms);
    azzero(fy, sys->natoms);
    azzero(fz, sys->natoms);

    /* precompute some constants */
    ssigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
    ssigma12 = ssigma6 * sigma * sigma * sigma * sigma * sigma * sigma;
    c12 = 4.0 * epsilon * ssigma12;
    c6 = 4.0 * epsilon * ssigma6;
    rcsq = rcut * rcut;

    double Box = sys->box;
    double halfbox = 0.5 * sys->box;

    for (int i = 0; i < ( sys->natoms) - 1; i += sys->nthreads) {
        int ii = i + tid;
        double rx, ry, rz;

        if (ii >= (sys->natoms - 1)) break;
        rx1 = sys->rx[ii];
        ry1 = sys->ry[ii];
        rz1 = sys->rz[ii];

        for (int j = ii + 1; j < (sys->natoms); ++j) {
            /* get distance between particle ii and j */
            rx = pbc(rx1 - sys->rx[j], halfbox, Box);
            ry = pbc(ry1 - sys->ry[j], halfbox, Box);
            rz = pbc(rz1 - sys->rz[j], halfbox, Box);
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

    #ifdef _OPENMP
    #pragma omp barrier
    #endif
    int i = (sys->natoms/sys->nthreads) + 1;
    int fromidx = tid * i;
    int toidx = fromidx + i;
    if (toidx > sys->natoms) toidx = sys->natoms;
    for (int i = 1; i < sys->nthreads; ++i){
      int offs = i*sys->natoms;
      for(int k = fromidx; k < toidx; ++k){
        sys->cx[k] += sys->cx[offs + k];
        sys->cy[k] += sys->cy[offs + k];
        sys->cz[k] += sys->cz[offs + k];
      }
    }
    
    for (int i = 0; i < sys->natoms; i++) sys->fx[i] = sys->cx[i];
    for (int i = 0; i < sys->natoms; i++) sys->fy[i] = sys->cy[i];
    for (int i = 0; i < sys->natoms; i++) sys->fz[i] = sys->cz[i];

    }

    sys->epot = epot;

}
