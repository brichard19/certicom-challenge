#ifndef _P131_CUH
#define _P131_CUH

#include "math_common.cuh"
#include "shared_types.h"

__constant__ uint131_t _p131_p = {{0x194c43186b3abc0b, 0x8e1d43f293469e33, 0x4}};
__constant__ uint131_t _p131_k = {{0xe0587d72985b105d, 0xf1fd54b0309e1ab9, 0x7cfd70cf}};
__constant__ uint131_t _p131_one = {{0x6e7743da32b6d0c7, 0x88c614d64c1a8f0b, 0x0}};
__constant__ uint131_t _p131_a = {{0xe7f7f250cee8709a, 0xacd15fe1a8ec1522, 0x0}};
__constant__ uint131_t _p131_b = {{0xc85087e5ab4eca9e, 0xde124657d7ba5851, 0x2}};

template <> struct Curve<131> {
  __device__ static uint131_t p() { return _p131_p; };
  __device__ static uint131_t one() { return _p131_one; };
  __device__ static uint131_t a() { return _p131_a; };
  __device__ static uint131_t b() { return _p131_b; };

  __device__ static uint131_t sub(uint131_t x, uint131_t y);
  __device__ static uint131_t add(uint131_t x, uint131_t y);
  __device__ static uint131_t mul(uint131_t x, uint131_t y);
};

__device__ uint131_t Curve<131>::sub(uint131_t x, uint131_t y)
{
  uint131_t z;

  int borrow = 0;

  uint64_t diff = x.w.v0 - y.w.v0;
  z.w.v0 = diff;
  borrow = diff > x.w.v0 ? 1 : 0;

  diff = x.w.v1 - y.w.v1 - borrow;
  z.w.v1 = diff;
  borrow = diff > x.w.v1 ? 1 : 0;

  // High word is only 3 bits
  diff = x.w.v2 - y.w.v2 - borrow;
  z.w.v2 = diff;
  borrow = diff & 0x08;

  // Went below zero. Need to add P.
  if(borrow) {
    int carry = 0;
    uint64_t sum = z.w.v0 + p().w.v0;
    z.w.v0 = sum;
    carry = sum < p().w.v0 ? 1 : 0;

    sum = z.w.v1 + p().w.v1 + carry;
    z.w.v1 = sum;
    carry = sum < p().w.v1 ? 1 : 0;

    z.w.v2 = z.w.v2 + p().w.v2 + carry;
  }

  return z;
}

__device__ uint131_t Curve<131>::add(uint131_t x, uint131_t y)
{
  uint131_t z = add_raw(x, y);

  // Reduce mod P
  if(is_less_than(p(), z)) {
    z = sub_raw(z, p());
  }
  return z;
}

// CIOS Montgomery multiplication for P131
// n=5 limbs of 32 bits, R=2^160, mp = -p^{-1} mod 2^32
// p.v[4] = 4 is hardcoded as a left-shift to save one multiply per iteration.
__device__ uint131_t Curve<131>::mul(uint131_t a, uint131_t b)
{
  const uint32_t mp = _p131_k.v[0]; // = 0x985b105d, -p^{-1} mod 2^32
  const uint32_t p0 = _p131_p.v[0];
  const uint32_t p1 = _p131_p.v[1];
  const uint32_t p2 = _p131_p.v[2];
  const uint32_t p3 = _p131_p.v[3];
  // p.v[4] = 4, use shift instead of multiply

  uint64_t t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0;

  for(int i = 0; i < 5; i++) {
    const uint32_t bi = b.v[i];
    uint64_t C, prod;

    // Multiply-accumulate: t += a * b[i]
    prod = (uint64_t)a.v[0] * bi + t0;
    t0 = prod & 0xffffffffULL;
    C = prod >> 32;
    prod = (uint64_t)a.v[1] * bi + t1 + C;
    t1 = prod & 0xffffffffULL;
    C = prod >> 32;
    prod = (uint64_t)a.v[2] * bi + t2 + C;
    t2 = prod & 0xffffffffULL;
    C = prod >> 32;
    prod = (uint64_t)a.v[3] * bi + t3 + C;
    t3 = prod & 0xffffffffULL;
    C = prod >> 32;
    // a.v[4] <= 7 (3 bits), so t4 stays well within uint64_t
    t4 = (uint64_t)a.v[4] * bi + t4 + C;

    // Montgomery reduction: m = t0 * mp mod 2^32, then t += m*p, shift right 32
    const uint32_t m = (uint32_t)t0 * mp;

    prod = (uint64_t)m * p0 + t0;
    C = prod >> 32;
    prod = (uint64_t)m * p1 + t1 + C;
    t0 = prod & 0xffffffffULL;
    C = prod >> 32;
    prod = (uint64_t)m * p2 + t2 + C;
    t1 = prod & 0xffffffffULL;
    C = prod >> 32;
    prod = (uint64_t)m * p3 + t3 + C;
    t2 = prod & 0xffffffffULL;
    C = prod >> 32;
    // p.v[4] = 4, so m * p.v[4] = m << 2
    prod = ((uint64_t)m << 2) + t4 + C;
    t3 = prod & 0xffffffffULL;
    t4 = prod >> 32;
  }

  uint131_t result;
  result.v[0] = (uint32_t)t0;
  result.v[1] = (uint32_t)t1;
  result.v[2] = (uint32_t)t2;
  result.v[3] = (uint32_t)t3;
  result.v[4] = (uint32_t)t4;

  if(is_less_than(_p131_p, result)) {
    result = sub_raw(result, _p131_p);
  }

  return result;
}

#endif