#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "ljmd.h"

extern const double kboltz;
extern const double mvsq2e;

/* helper function: get current time in seconds since epoch */
double wallclock()
{
        struct timeval t;
        gettimeofday(&t,0);
        return ((double) t.tv_sec) + 1.0e-6*((double) t.tv_usec);
}


/* compute kinetic energy */
void ekin(mdsys_t *sys)
{
    int i;
    double natoms = sys->natoms;
    double mass = sys->mass;
    double *vx = sys->vx;
    double *vy = sys->vy;
    double *vz = sys->vz;
    double constant = 0.5 * mvsq2e * mass;
    double constant2 = 2.0 / (3.0*natoms-3.0) * 1.0 / kboltz;
    sys->ekin=0.0;
    for (i=0; i<natoms; ++i) {
        sys->ekin += constant * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    }
    sys->temp = sys->ekin * constant2;
}

allocate_mem(mdsys_t *sys){

    sys->rx=(double *)malloc(sys->natoms*sizeof(double));
    sys->ry=(double *)malloc(sys->natoms*sizeof(double));
    sys->rz=(double *)malloc(sys->natoms*sizeof(double));
    sys->vx=(double *)malloc(sys->natoms*sizeof(double));
    sys->vy=(double *)malloc(sys->natoms*sizeof(double));
    sys->vz=(double *)malloc(sys->natoms*sizeof(double));
    sys->fx=(double *)malloc(sys->natoms*sizeof(double));
    sys->fy=(double *)malloc(sys->natoms*sizeof(double));
    sys->fz=(double *)malloc(sys->natoms*sizeof(double));

    //allocate memory for auxiliary
    sys->cx=(double *)malloc(sys->natoms*sizeof(double));
    sys->cy=(double *)malloc(sys->natoms*sizeof(double));
    sys->cz=(double *)malloc(sys->natoms*sizeof(double));

}
