#ifndef _SHARED_TYPES_H
#define _SHARED_TYPES_H

#include <stdint.h>

// Types that are used in both host and device code

union uint131_t {
  struct {
    uint64_t v0;
    uint64_t v1;
    uint32_t v2;
  } w;
  uint32_t v[5];
};

union uint160_t {
  struct {
    uint64_t v0;
    uint64_t v1;
    uint32_t v2;
  } w;
};

struct uint262_t {
  uint64_t v[5];
};

typedef unsigned __int128 uint128_t;

#endif