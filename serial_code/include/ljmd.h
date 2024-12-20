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
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
};
typedef struct _mdsys mdsys_t;

extern double wallclock();

/* helper function: zero out an array */
extern void azzero(double *d, const int n);

/* helper function: apply minimum image convention */
extern double pbc(double x, const double boxby2);

extern void ekin(mdsys_t *sys);

extern void force(mdsys_t *sys);

extern void velverlet_step1(mdsys_t *sys);

extern void velverlet_step2(mdsys_t *sys);

extern void output(mdsys_t *sys, FILE *erg, FILE *traj);

extern int get_a_line(FILE *fp, char *buf);

extern void reading(char *ll, mdsys_t *sys, char *restfile, char *trajfile, char *ergfile, int *nprint);

extern void allocation(mdsys_t *sys);

extern void cleanup(mdsys_t *sys);

#endif
