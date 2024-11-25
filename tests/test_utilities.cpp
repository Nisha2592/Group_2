#include "gtest/gtest.h"
#include "ljmd.h"
#include <cmath>

TEST(TestAzzero, doubles) {
  double *buf = new double[10];
  for (int i = 0; i < 10; i++) buf[i] = i +1;
  ASSERT_DOUBLE_EQ(buf[1], 2.0);
  ASSERT_DOUBLE_EQ(buf[5], 6.0);
  ASSERT_DOUBLE_EQ(buf[9], 10.0);

  azzero(buf, 10);
  ASSERT_DOUBLE_EQ(buf[1], 0.0);
  ASSERT_DOUBLE_EQ(buf[5], 0.0);
  ASSERT_DOUBLE_EQ(buf[9], 0.0);
}

class PBCTest:
  public::testing::Test {
    protected:
    mdsys_t * sys;
    void SetUp() {
      sys = new mdsys_t;
      sys->natoms =2;
      sys->rx = new double[2];
      sys->ry = new double[2];
      sys->rz = new double[2];
      sys->rx[0] = -6.0;
      sys->rx[1] =  6.0;
      sys->ry[0] = -6.0;
      sys->ry[1] =  6.0;
      sys->rz[0] = -6.0;
      sys->rz[1] =  6.0;
      sys->box = 10;
  }

  void TearDown() {
    delete[] sys->rx;
    delete[] sys->ry;
    delete[] sys->rz;

    delete sys;
  }
};

TEST_F(PBCTest, doubles) {
   int i = 0;
   int j = 1;
   ASSERT_NE(sys,nullptr);
   double rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
   double ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
   double rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
   ASSERT_DOUBLE_EQ(rx*rx + ry*ry + rz*rz, 12);
}
