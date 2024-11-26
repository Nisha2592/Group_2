#ifndef C_HEADERS

#define C_HEADERS 

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


/* generic file- or pathname buffer length */
#define BLEN 200

extern const double kboltz;
extern const double mvsq2e;

/* structure to hold the complete information
 * about the MD system */
struct _mdsys{
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    int natoms,nfi,nsteps;
};
typedef struct _mdsys mdsys_t;

extern double wallclock();

/* helper function: zero out an array */
static __attribute__((always_inline))
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* helper function: apply minimum image convention */
static __attribute__((always_inline,pure))
double pbc(double x, const double boxby2, const double Box)
{
    while (x >  boxby2) x -= Box;
    while (x < -boxby2) x += Box;
    return x;
}

extern void ekin(mdsys_t *sys);

extern void force(mdsys_t *sys);

extern void velverlet_step1(mdsys_t *sys);

extern void velverlet_step2(mdsys_t *sys);

extern void output(mdsys_t *sys, FILE *erg, FILE *traj);

extern int get_a_line(FILE *fp, char *buf);

extern void cleanup(mdsys_t *sys);

#endif
