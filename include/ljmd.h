#ifndef C_HEADERS

#define C_HEADERS 

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>


/* generic file- or pathname buffer length */
#define BLEN 200

extern const double kboltz;
extern const double mvsq2e;

/* structure to hold the complete information
 * about the MD system */
struct _mdsys{ 
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    double *cx, *cy, *cz;
    
};
typedef struct _mdsys mdsys_t;

extern double wallclock();

extern void azzero(double *d, const int n);

extern double pbc(double x, const double boxby2);

extern void ekin(mdsys_t *sys);

extern void force(mdsys_t *sys, int rank, int size, MPI_Comm comm);

extern void velverlet_step1(mdsys_t *sys);

extern void velverlet_step2(mdsys_t *sys);

extern void output(mdsys_t *sys, FILE *erg, FILE *traj);

extern int get_a_line(FILE *fp, char *buf);

extern void allocate_mem(mdsys_t *sys);

extern void cleanup(mdsys_t *sys);


#endif
