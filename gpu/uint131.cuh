#ifndef _UINT131_CUH
#define _UINT131_CUH

#include <stdint.h>

// Represented a 131-bit integer
struct uint131_t {
  uint64_t v[3];
};

#endif