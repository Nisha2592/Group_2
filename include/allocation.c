#include <stdio.h>
#include <ctype.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "ljmd.h"

void allocation(mdsys_t *sys)
{
    sys->rx=(double *)malloc(sys->natoms*sizeof(double));
    sys->ry=(double *)malloc(sys->natoms*sizeof(double));
    sys->rz=(double *)malloc(sys->natoms*sizeof(double));

    sys->vx=(double *)malloc(sys->natoms*sizeof(double));
    sys->vy=(double *)malloc(sys->natoms*sizeof(double));
    sys->vz=(double *)malloc(sys->natoms*sizeof(double));

    sys->fx=(double *)malloc(sys->natoms*sizeof(double));
    sys->fy=(double *)malloc(sys->natoms*sizeof(double));
    sys->fz=(double *)malloc(sys->natoms*sizeof(double));

    #if defined(_OPENMP)
    #pragma omp parallel
    {
    if(0 == omp_get_thread_num()){
      sys->nthreads = omp_get_num_threads();
    }
    }
    #else 
    sys->nthreads = 1;
    #endif
    
    //allocate memory for enlarged auxilliary buffers 
    sys->cx=(double *)malloc(sys->nthreads*sys->natoms*sizeof(double));
    sys->cy=(double *)malloc(sys->nthreads*sys->natoms*sizeof(double));
    sys->cz=(double *)malloc(sys->nthreads*sys->natoms*sizeof(double));
}