#include "gtest/gtest.h"
#include "../include/ljmd.h"
#include "../include/utilities.c"
#include "../include/force.c"

class VerletTest:
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
      sys->vx[0] = 0.0;
      sys->vx[1] = 0.0;
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
    void Free() {
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
  TEST_F(VerletTest, step1) {
    ASSERT_NE(sys, nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
    ASSERT_DOUBLE_EQ(sys->vx[0], 0.0);
    velverlet_step1(sys);
    ASSERT_DOUBLE_EQ(sys->rx[0], 0.5*mvsq2e - 1);
    ASSERT_DOUBLE_EQ(sys->vx[0], 0.5);
  }

  TEST_F(VerletTest, step2) {
    ASSERT_NE(sys,nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
    ASSERT_DOUBLE_EQ(sys->vx[0], 0.0);
    velverlet_step2(sys);
    ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
    ASSERT_DOUBLE_EQ(sys->vx[0], 0.5);
  }
 
