#ifndef _P89_CUH
#define _P89_CUH

#include "math_common.cuh"
#include "shared_types.h"

__constant__ uint131_t _p89_p = {{0x903f1643908ba955, 0x158685c, 0x0}};
__constant__ uint131_t _p89_k = {{0x40a0dbc30d18f403, 0x5ed01b73e78a4d11, 0xf15c3823}};
__constant__ uint131_t _p89_r2 = {{0xd96623186da369dd, 0x13e5217, 0x0}};
__constant__ uint131_t _p89_one = {{0x9e896dc839826cf6, 0x78c8ca, 0x0}};
__constant__ uint131_t _p89_a = {{0xb1ec24706b2573c1, 0x593048, 0x0}};
__constant__ uint131_t _p89_b = {{0x37717b7e1c0a25af, 0xc58c6, 0x0}};

template <> struct Curve<89> {
  __device__ static uint131_t p() { return _p89_p; };
  __device__ static uint131_t k() { return _p89_k; };
  __device__ static uint131_t r2() { return _p89_r2; };
  __device__ static uint131_t one() { return _p89_one; };
  __device__ static uint131_t a() { return _p89_a; };
  __device__ static uint131_t b() { return _p89_b; };

  __device__ static uint131_t sub(uint131_t x, uint131_t y);
  __device__ static uint131_t add(uint131_t x, uint131_t y);
  __device__ static uint131_t fused_mul(uint131_t x, uint131_t y);
};

__device__ uint131_t Curve<89>::sub(uint131_t x, uint131_t y)
{
  uint131_t z = {{0}};

  int borrow = 0;

  uint64_t diff = x.w.v0 - y.w.v0;
  z.w.v0 = diff;
  borrow = diff > x.w.v0 ? 1 : 0;

  diff = x.w.v1 - y.w.v1 - borrow;
  z.w.v1 = diff;
  borrow = diff & (1 << 25);

  //   // High word is only 3 bits
  //   diff = x.w.v2 - y.w.v2 - borrow;
  //   z.w.v2 = diff;
  //   borrow = diff & 0x08;

  // Went below zero. Need to add P.
  if(borrow) {
    int carry = 0;
    uint64_t sum = z.w.v0 + p().w.v0;
    z.w.v0 = sum;
    carry = sum < p().w.v0 ? 1 : 0;

    sum = z.w.v1 + p().w.v1 + carry;
    z.w.v1 = sum;
    // carry = sum < _p89_p.w.v1 ? 1 : 0;

    // z.w.v2 = z.w.v2 + _p89_p.w.v2 + carry;
  }

  return z;
}

__device__ uint131_t Curve<89>::add(uint131_t x, uint131_t y)
{
  uint131_t z = add_raw(x, y);

  // Reduce mod P
  if(is_less_than(p(), z)) {
    z = sub_raw(z, p());
  }
  return z;
}

// CIOS Montgomery multiplication for P89
// n=5 limbs of 32 bits, R=2^160, mp = -p^{-1} mod 2^32
// p has only 3 significant limbs: p.v[3]=p.v[4]=0
// b has only 3 significant limbs: b.v[3]=b.v[4]=0 (values < 2^89)
// Outer loop i=3,4: multiply step is zero; still need reduction passes.
__device__ uint131_t Curve<89>::fused_mul(uint131_t a, uint131_t b)
{
  const uint32_t mp = _p89_k.v[0]; // = 0x0d18f403, -p^{-1} mod 2^32
  const uint32_t p0 = _p89_p.v[0]; // = 0x908ba955
  const uint32_t p1 = _p89_p.v[1]; // = 0x903f1643
  const uint32_t p2 = _p89_p.v[2]; // = 0x0158685c
  // p.v[3] = p.v[4] = 0

  uint64_t t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0;

  // i = 0, 1, 2: full multiply-accumulate + reduction
  for(int i = 0; i < 3; i++) {
    const uint32_t bi = b.v[i];
    uint64_t C, prod;

    // Multiply-accumulate: t += a[0..2] * b[i]  (a.v[3]=a.v[4]=0)
    prod = (uint64_t)a.v[0] * bi + t0;
    t0 = prod & 0xffffffffULL;
    C = prod >> 32;
    prod = (uint64_t)a.v[1] * bi + t1 + C;
    t1 = prod & 0xffffffffULL;
    C = prod >> 32;
    prod = (uint64_t)a.v[2] * bi + t2 + C;
    t2 = prod & 0xffffffffULL;
    C = prod >> 32;
    t3 = t3 + C;

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
    // p.v[3] = 0
    t2 = t3 + C;
    t3 = t4;
    t4 = 0;
  }

  // i = 3, 4: b[i] = 0 so multiply step is zero; only reduction needed
  for(int i = 3; i < 5; i++) {
    uint64_t C, prod;

    const uint32_t m = (uint32_t)t0 * mp;

    prod = (uint64_t)m * p0 + t0;
    C = prod >> 32;
    prod = (uint64_t)m * p1 + t1 + C;
    t0 = prod & 0xffffffffULL;
    C = prod >> 32;
    prod = (uint64_t)m * p2 + t2 + C;
    t1 = prod & 0xffffffffULL;
    C = prod >> 32;
    // p.v[3] = 0
    t2 = t3 + C;
    t3 = t4;
    t4 = 0;
  }

  uint131_t result;
  result.v[0] = (uint32_t)t0;
  result.v[1] = (uint32_t)t1;
  result.v[2] = (uint32_t)t2;
  result.v[3] = (uint32_t)t3;
  result.v[4] = (uint32_t)t4;

  if(is_less_than(_p89_p, result)) {
    result = sub_raw(result, _p89_p);
  }

  return result;
}

#endif