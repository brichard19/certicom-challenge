#ifndef _UINT131_CUH
#define _UINT131_CUH

#include <stdint.h>

// Represented a 131-bit integer
union uint131_t {
  struct {
    uint64_t v0;
    uint64_t v1;
    uint32_t v2;
  }w;
};

union uint160_t {
  struct {
    uint64_t v0;
    uint64_t v1;
    uint32_t v2;
  }w;
};

struct uint262_t {
  uint64_t v[5];
};

#endif