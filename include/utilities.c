#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "ljmd.h"


//const double kboltz=0.0019872067;
//const double mvsq2e=2390.05736153349;
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
