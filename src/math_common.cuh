#ifndef _MATH_COMMON_CUH
#define _MATH_COMMON_CUH

#include "shared_types.h"

template<int CURVE> struct Curve {};

// Memory layout for array of 131-bit integers
//
// Each 131-bit integer is 17 bytes. Stored in memory as 1 x 16-byte element and 1 x 1-byte
//
// For coalesed memory access, the 2 parts are stored in separate contiguous regions of memmory
//
// Accessing array element at index 'idx' in an array size 'n':
//
// Byte offset for 16-byte portion:
//
// idx * 16
//
// Byte offset for 1-byte portion:
//
// n * 16 + idx
//
typedef unsigned __int128 uint128_t;


__device__ uint131_t load_uint131(const void* p, int idx, int n)
{
  const uint128_t* p128 = (const uint128_t*)p;
  const uint8_t* p8 = (const uint8_t*)p;

  uint128_t u128 = p128[idx];
  
  uint131_t x;
  x.w.v0 = (uint64_t)u128;
  x.w.v1 = (uint64_t)(u128 >> 64);
  x.w.v2 = (uint32_t)p8[n * sizeof(uint128_t) + idx];

  return x;
}

__device__ void store_uint131(void* p, int idx, int n, uint131_t x)
{
  uint128_t* p128 = (uint128_t*)p;
  uint8_t* p8 = (uint8_t*)p;

  uint128_t u128 = ((uint128_t)x.w.v1 << 64) | x.w.v0;
  p128[idx] = u128;

  p8[n * sizeof(uint128_t) + idx] = (uint8_t)x.w.v2;
}


__device__ uint131_t sub_raw(const uint131_t& x, const uint131_t& y)
{
  uint131_t z;
  int borrow = 0;

  uint128_t diff = (uint128_t)x.w.v0 - y.w.v0 - borrow;
  z.w.v0 = (uint64_t)diff;
  borrow = (int)(diff >> 64) & 1;
  
  diff = (uint128_t)x.w.v1 - y.w.v1 - borrow;
  z.w.v1 = (uint64_t)diff;
  borrow = (int)(diff >> 64) & 1;

  z.w.v2 = x.w.v2 - y.w.v2 - borrow;

  return z;
}

__device__ uint131_t add_raw(const uint131_t& x, const uint131_t& y, int carry_in = 0)
{

  uint131_t z;
  int carry = carry_in;
  
  uint128_t sum = (uint128_t)x.w.v0 + y.w.v0 + carry;
  z.w.v0 = (uint64_t)sum;
  carry = (uint64_t)(sum >> 64);
  
  sum = (uint128_t)x.w.v1 + y.w.v1 + carry;
  z.w.v1 = (uint64_t)sum;
  carry = (uint64_t)(sum >> 64);
 
  z.w.v2 = x.w.v2 + y.w.v2 + carry;

  return z;
}

__device__ int is_less_than(uint131_t x, uint131_t y)
{
  uint131_t diff = sub_raw(x, y);

  return (diff.w.v2 >> 31) & 1;
}

__device__ bool equal(const uint131_t& x, const uint131_t& y)
{
  return x.w.v0 == y.w.v0 && x.w.v1 == y.w.v1 && x.w.v2 == y.w.v2;
}



__device__ uint160_t mul_mod_160(const uint160_t& a, const uint131_t& b)
{
  uint160_t tmp;
  uint64_t high64 = 0;
  uint32_t high32 = 0;

  // a0 * b0
  uint128_t t = (uint128_t)a.w.v0 * b.w.v0;
  tmp.w.v0 = (uint64_t)t;
  high64 = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a.w.v0 * b.w.v1 + high64;
  tmp.w.v1 = (uint64_t)t;
  high32 = (uint32_t)(t >> 64);

  // a0 * b2
  uint32_t t32 = (uint32_t)a.w.v0 * b.w.v2 + high32;
  tmp.w.v2 = (uint32_t)t32;

  // a1 * b0
  t = (uint128_t)a.w.v1 * b.w.v0 + tmp.w.v1;
  tmp.w.v1 = (uint64_t)t;
  high32 = (uint32_t)(t >> 64);

  // a0 * b1 
  t32 = (uint32_t)a.w.v1 * b.w.v1 + tmp.w.v2 + high32;
  tmp.w.v2 = t32;

  // a2 * b0
  t32 = (uint32_t)a.w.v2 * b.w.v0 + tmp.w.v2;
  tmp.w.v2 = t32;

  return tmp;
}

#endif