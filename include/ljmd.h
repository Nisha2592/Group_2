#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
};
typedef struct _mdsys mdsys_t;

static double wallclock();

static void azzero(double *d, const int n);

static double pbc(double x, const double boxby2);

static void ekin(mdsys_t *sys);

static void force(mdsys_t *sys);

static void velverlet(mdsys_t *sys);

static void output(mdsys_t *sys, FILE *erg, FILE *traj);