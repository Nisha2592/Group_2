#include <math.h>
#include <mpi.h>
#include <sys/time.h>
#include "ljmd.h"

const double kboltz=0.0019872067;
const double mvsq2e=2390.05736153349;

void velverlet_step1(mdsys_t *sys) {
    int i;
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        sys->rx[i] += sys->dt * sys->vx[i];
        sys->ry[i] += sys->dt * sys->vy[i];
        sys->rz[i] += sys->dt * sys->vz[i];
    }

}

void velverlet_step2(mdsys_t *sys) {
    /* compute forces and potential energy */
    //force(sys);
    int i;
    for (i = 0; i < sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }
}

void force(mdsys_t *sys, int rank, int size, MPI_Comm comm) {
    double epot_local = 0.0;  // Local potential energy for each rank

    // Initialize forces to zero
    double epot =0.0;
    azzero(sys->cx, sys->natoms);
    azzero(sys->cy, sys->natoms);
    azzero(sys->cz, sys->natoms);

    // Broadcast particle positions to all ranks
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Get rank and size
//    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calculate chunk size for each rank
    int chunk_size = sys->natoms / size;
    int remainder = sys->natoms % size;

    // Calculate start and end indices for each rank
    int start_idx = rank * chunk_size + (rank < remainder ? rank : remainder);
    int end_idx = start_idx + chunk_size + (rank < remainder ? 1 : 0);

    // Precompute constants for force and potential energy calculations
    double rcsq = sys->rcut * sys->rcut;
    double c6 = sys->sigma * sys->sigma * sys->sigma * sys->sigma * sys->sigma * sys->sigma;
    double c12 = 4.0 * sys->epsilon * c6 * c6;
    c6 *= 4.0 * sys->epsilon;

    // Loop over assigned atoms (Newton's third law applied here)
    for (int i = start_idx; i < end_idx; ++i) {
        double rx_i = sys->rx[i];
        double ry_i = sys->ry[i];
        double rz_i = sys->rz[i];

        for (int j = i + 1; j < sys->natoms; ++j) {
            // Calculate distance between atoms i and j
            double rx = pbc(rx_i - sys->rx[j], 0.5 * sys->box, sys->box);
            double ry = pbc(ry_i - sys->ry[j], 0.5 * sys->box, sys->box);
            double rz = pbc(rz_i - sys->rz[j], 0.5 * sys->box, sys->box);
            double rsq = rx * rx + ry * ry + rz * rz;

            // Apply cutoff distance
            if (rsq < rcsq) {
                double rinv = 1.0 / rsq;
                double rinv6 = rinv * rinv * rinv;
                double rinv12 = rinv6 * rinv6;
                double ffac = (12.0 * c12 * rinv12 - 6.0 * c6 * rinv6) * rinv;
                epot_local += (c12 * rinv12 - c6 * rinv6);

                // Update forces
                sys->cx[i] += rx * ffac;
                sys->cy[i] += ry * ffac;
                sys->cz[i] += rz * ffac;

                sys->cx[j] -= rx * ffac;  // Newton's third law
                sys->cy[j] -= ry * ffac;  // Newton's third law
                sys->cz[j] -= rz * ffac;  // Newton's third law
            }
        }
    }

    // Reduce forces from all ranks into global arrays
    MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Reduce potential energy across ranks
    MPI_Reduce(&epot_local, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

