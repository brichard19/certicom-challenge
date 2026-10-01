#ifndef _P109_CUH
#define _P109_CUH

#include "math_common.cuh"
#include "shared_types.h"

__constant__ uint131_t _p109_p = {{0x0b5b521e6d9fb599, 0x00001bd579792b38, 0x0}};
__constant__ uint131_t _p109_k = {{0x95f4da11710f3157, 0x2a9b9beb62e06384, 0x7bc599a4}};
__constant__ uint131_t _p109_one = {{0xb2ec742432d0a6bf, 0x00000accf5a4e730, 0x0}};
__constant__ uint131_t _p109_a = {{0x0bc15e225a5ee910, 0x00000dc462a14129, 0x0}};
__constant__ uint131_t _p109_b = {{0xb52890719eff9942, 0x0000052e3b8cfeee, 0x0}};

template <> struct Curve<CURVE_ID_ECP109> {
  __device__ static uint131_t p() { return _p109_p; };
  __device__ static uint131_t one() { return _p109_one; };
  __device__ static uint131_t a() { return _p109_a; };
  __device__ static uint131_t b() { return _p109_b; };

  __device__ static uint131_t sub(uint131_t x, uint131_t y);
  __device__ static uint131_t add(uint131_t x, uint131_t y);
  __device__ static uint131_t mul(uint131_t x, uint131_t y);
};

__device__ uint131_t Curve<CURVE_ID_ECP109>::sub(uint131_t x, uint131_t y)
{
  uint131_t z = {{0}};

  uint64_t diff = x.w.v0 - y.w.v0;
  uint64_t borrow = diff > x.w.v0;
  z.w.v0 = diff;

  diff = x.w.v1 - y.w.v1 - borrow;
  z.w.v1 = diff;

  // Valid field elements use 45 bits of the high word. An underflow sets bit 45.
  if(diff & (uint64_t(1) << 45)) {
    uint64_t sum = z.w.v0 + p().w.v0;
    uint64_t carry = sum < p().w.v0;
    z.w.v0 = sum;
    z.w.v1 += p().w.v1 + carry;
  }

  return z;
}

__device__ uint131_t Curve<CURVE_ID_ECP109>::add(uint131_t x, uint131_t y)
{
  uint131_t z = add_raw(x, y);

  if(!is_less_than(z, p())) {
    z = sub_raw(z, p());
  }
  return z;
}

// CIOS Montgomery multiplication for P109.
// R=2^160 with four significant 32-bit limbs; limb 4 is zero.
__device__ uint131_t Curve<CURVE_ID_ECP109>::mul(uint131_t a, uint131_t b)
{
  const uint32_t mp = _p109_k.v[0]; // -p^{-1} mod 2^32
  const uint32_t p0 = _p109_p.v[0];
  const uint32_t p1 = _p109_p.v[1];
  const uint32_t p2 = _p109_p.v[2];
  const uint32_t p3 = _p109_p.v[3];

  uint64_t t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0;

  // Four nonzero limbs of b require multiply-accumulate and reduction.
  for(int i = 0; i < 4; i++) {
    const uint32_t bi = b.v[i];
    uint64_t carry, product;

    product = (uint64_t)a.v[0] * bi + t0;
    t0 = (uint32_t)product;
    carry = product >> 32;
    product = (uint64_t)a.v[1] * bi + t1 + carry;
    t1 = (uint32_t)product;
    carry = product >> 32;
    product = (uint64_t)a.v[2] * bi + t2 + carry;
    t2 = (uint32_t)product;
    carry = product >> 32;
    product = (uint64_t)a.v[3] * bi + t3 + carry;
    t3 = (uint32_t)product;
    t4 += product >> 32;

    const uint32_t m = (uint32_t)t0 * mp;
    product = (uint64_t)m * p0 + t0;
    carry = product >> 32;
    product = (uint64_t)m * p1 + t1 + carry;
    t0 = (uint32_t)product;
    carry = product >> 32;
    product = (uint64_t)m * p2 + t2 + carry;
    t1 = (uint32_t)product;
    carry = product >> 32;
    product = (uint64_t)m * p3 + t3 + carry;
    t2 = (uint32_t)product;
    carry = product >> 32;
    t3 = t4 + carry;
    t4 = 0;
  }

  // b.v[4] is zero, leaving one final Montgomery reduction pass.
  {
    uint64_t carry, product;
    const uint32_t m = (uint32_t)t0 * mp;

    product = (uint64_t)m * p0 + t0;
    carry = product >> 32;
    product = (uint64_t)m * p1 + t1 + carry;
    t0 = (uint32_t)product;
    carry = product >> 32;
    product = (uint64_t)m * p2 + t2 + carry;
    t1 = (uint32_t)product;
    carry = product >> 32;
    product = (uint64_t)m * p3 + t3 + carry;
    t2 = (uint32_t)product;
    carry = product >> 32;
    t3 = t4 + carry;
    t4 = 0;
  }

  uint131_t result = {{0}};
  result.v[0] = (uint32_t)t0;
  result.v[1] = (uint32_t)t1;
  result.v[2] = (uint32_t)t2;
  result.v[3] = (uint32_t)t3;

  if(!is_less_than(result, _p109_p)) {
    result = sub_raw(result, _p109_p);
  }
  return result;
}

#endif
