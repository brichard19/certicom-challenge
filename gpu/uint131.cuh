#ifndef _UINT131_CUH
#define _UINT131_CUH

#include <stdint.h>

// Represented a 131-bit integer
union uint131_t {
  uint64_t v[3];
  struct {
    uint64_t v0;
    uint64_t v1;
    uint32_t v2;
    uint32_t padding;
  }w;
};

union uint160_t {
  uint64_t v[3];
  struct {
    uint64_t v0;
    uint64_t v1;
    uint32_t v2;
    uint32_t padding;
  }w;
};

struct uint262_t {
  uint64_t v[5];
};

#endif