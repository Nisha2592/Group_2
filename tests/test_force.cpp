#include "gtest/gtest.h"
#include "ljmd.h"
#include "comp.h"


class forceTest:
  public::testing::Test {
    protected:
    mdsys_t *sys;
    void SetUp() {
      sys = new mdsys_t;
      sys->natoms = 2;
      sys->mass = 1.0;
      sys->ekin = 0.0;
      sys->temp = 0.0;
      sys->box = 10.0;
      sys->rcut = 2.5;
      sys->epsilon = 1.0;
      sys->sigma = 1.0;
      sys->dt = mvsq2e;
      sys->rx = new double[2];
      sys->vx = new double[2];
      sys->fx = new double[2];
      sys->ry = new double[2];
      sys->vy = new double[2];
      sys->fy = new double[2];
      sys->rz = new double[2];
      sys->vz = new double[2];
      sys->fz = new double[2];
      sys->rx[0] = -1.0;
      sys->rx[1] = 1.0;
      sys->ry[0] = -1.0;
      sys->ry[1] = 1.0;
      sys->rz[0] = -1.0;
      sys->rz[1] = 0.0;
      sys->vx[0] = 1.0;
      sys->vx[1] = 1.0;
      sys->vy[0] = 0.0;
      sys->vy[1] = 0.0;
      sys->vz[0] = 0.0;
      sys->vz[1] = 0.0;
      sys->fx[0] = 1.0;
      sys->fx[1] = 0.2;
      sys->fy[0] = 1.0;
      sys->fy[1] = 0.2;
      sys->fz[0] = 1.0;
      sys->fz[1] = 0.2;
    }
    void TearDown() {
      delete[] sys->rx;
      delete[] sys->vx;
      delete[] sys->fx;
      delete[] sys->ry;
      delete[] sys->vy;
      delete[] sys->fy;
      delete[] sys->rz;
      delete[] sys->vz;
      delete[] sys->fz;

      delete sys;
    }
  };

// Test for force calculation inside and outside the cutoff distance
TEST_F(forceTest, ForceCalculation) {
    // Setup particles inside the cutoff
    sys->rx[0] = 0.0;
    sys->ry[0] = 0.0;
    sys->rz[0] = 0.0;
    sys->rx[1] = 1.0; // Within the cutoff distance
    sys->ry[1] = 0.0;
    sys->rz[1] = 0.0;

    // Call the force calculation
    force(sys);

    // Check potential energy and force calculation inside the cutoff
    double expected_potential_energy = 4.0 * (pow(1.0 / 1.0, 12) - pow(1.0 / 1.0, 6));
    ASSERT_NEAR(sys->epot, expected_potential_energy, 1e-5);

    // Check if forces are non-zero (as they should be inside the cutoff)
    ASSERT_NEAR(sys->fx[0], -sys->fx[1], 1e-5);
    ASSERT_NEAR(sys->fy[0], -sys->fy[1], 1e-5);
    ASSERT_NEAR(sys->fz[0], -sys->fz[1], 1e-5);

    // Test particles outside the cutoff distance
    sys->rx[1] = 3.0; // Outside the cutoff distance

    // Reset forces
    azzero(sys->fx, sys->natoms);
    azzero(sys->fy, sys->natoms);
    azzero(sys->fz, sys->natoms);

    // Call the force calculation again
    force(sys);

    // Forces should be near zero if particles are outside the cutoff
    ASSERT_NEAR(sys->fx[0], 0.0, 1e-10);
    ASSERT_NEAR(sys->fy[0], 0.0, 1e-10);
    ASSERT_NEAR(sys->fz[0], 0.0, 1e-10);
}

// Test for kinetic energy calculation

 TEST_F(forceTest, KineticEnergyCalculation) {
   // initialize velocities
   sys->vx[0] = 1.0;
   sys->vy[0] = 0.0;
   sys->vz[0] = 0.0;
   sys->vx[1] = 1.0;
   sys->vy[1] = 0.0;
   sys->vz[1] = 0.0;
   ASSERT_NE(sys, nullptr);
   ASSERT_DOUBLE_EQ(sys->vx[0], 1.0);
   ASSERT_DOUBLE_EQ(sys->vx[1], 1.0);
   ekin(sys);
   ASSERT_DOUBLE_EQ(sys->ekin, mvsq2e);
   ASSERT_DOUBLE_EQ(sys->temp, 2.0*mvsq2e/(3.0*kboltz));
}

