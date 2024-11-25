#include "ljmd.h"

#ifndef COMPUTE_H
  #define COMPUTE_H
  #ifdef __cplusplus
  extern "C" {
  #endif
  void ekin(mdsys_t *sys);
  void force(mdsys_t *sys);
  #ifdef __cplusplus
  }
  #endif
#endif
