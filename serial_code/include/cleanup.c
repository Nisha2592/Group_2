// cleanup.c

#include <stdlib.h>
#include "ljmd.h"

void cleanup(mdsys_t *sys) {
    free(sys->rx);
    free(sys->ry);
    free(sys->rz);
    free(sys->vx);
    free(sys->vy);
    free(sys->vz);
    free(sys->fx);
    free(sys->fy);
    free(sys->fz);
}

