#ifndef _P131_CUH
#define _P131_CUH

#include "math_common.cuh"
#include "shared_types.h"

__constant__ uint131_t _p131_p = {{0x194c43186b3abc0b, 0x8e1d43f293469e33, 0x4}};
__constant__ uint131_t _p131_k = {{0xe0587d72985b105d, 0xf1fd54b0309e1ab9, 0x7cfd70cf}};
__constant__ uint131_t _p131_one = {{0x6e7743da32b6d0c7, 0x88c614d64c1a8f0b, 0x0}};
__constant__ uint131_t _p131_a = {{0xe7f7f250cee8709a, 0xacd15fe1a8ec1522, 0x0}};
__constant__ uint131_t _p131_b = {{0xc85087e5ab4eca9e, 0xde124657d7ba5851, 0x2}};

template <> struct Curve<CURVE_ID_ECP131> {
  __device__ static uint131_t p() { return _p131_p; };
  __device__ static uint131_t one() { return _p131_one; };
  __device__ static uint131_t a() { return _p131_a; };
  __device__ static uint131_t b() { return _p131_b; };

  __device__ static uint131_t sub(uint131_t x, uint131_t y);
  __device__ static uint131_t add(uint131_t x, uint131_t y);
  __device__ static uint131_t mul(uint131_t x, uint131_t y);
  __device__ static uint131_t square(uint131_t x);
};

__device__ uint131_t Curve<CURVE_ID_ECP131>::sub(uint131_t x, uint131_t y)
{
  uint131_t z;

#if defined(__HIP_DEVICE_COMPILE__) && (defined(__GFX11__) || defined(__GFX12__))
  uint32_t carry, p0, p1, p2, p3, p4;
  // Keep the borrow/carry chains together, and add P only to underflowing lanes.
  // Early-clobber outputs must not overlap inputs consumed later in the chain.
  asm("v_sub_co_u32 %[z0], %[cc], %[x0], %[y0]\n\t"
      "v_sub_co_ci_u32 %[z1], %[cc], %[x1], %[y1], %[cc]\n\t"
      "v_sub_co_ci_u32 %[z2], %[cc], %[x2], %[y2], %[cc]\n\t"
      "v_sub_co_ci_u32 %[z3], %[cc], %[x3], %[y3], %[cc]\n\t"
      "v_sub_co_ci_u32 %[z4], %[cc], %[x4], %[y4], %[cc]\n\t"
      "v_cndmask_b32 %[p0], 0, %[mod0], %[cc]\n\t"
      "v_cndmask_b32 %[p1], 0, %[mod1], %[cc]\n\t"
      "v_cndmask_b32 %[p2], 0, %[mod2], %[cc]\n\t"
      "v_cndmask_b32 %[p3], 0, %[mod3], %[cc]\n\t"
      "v_cndmask_b32 %[p4], 0, %[mod4], %[cc]\n\t"
      "v_add_co_u32 %[z0], %[cc], %[z0], %[p0]\n\t"
      "v_add_co_ci_u32 %[z1], %[cc], %[z1], %[p1], %[cc]\n\t"
      "v_add_co_ci_u32 %[z2], %[cc], %[z2], %[p2], %[cc]\n\t"
      "v_add_co_ci_u32 %[z3], %[cc], %[z3], %[p3], %[cc]\n\t"
      "v_add_co_ci_u32 %[z4], %[cc], %[z4], %[p4], %[cc]"
      : [z0] "=&v"(z.v[0]), [z1] "=&v"(z.v[1]), [z2] "=&v"(z.v[2]),
        [z3] "=&v"(z.v[3]), [z4] "=&v"(z.v[4]), [cc] "=&s"(carry),
        [p0] "=&v"(p0), [p1] "=&v"(p1), [p2] "=&v"(p2),
        [p3] "=&v"(p3), [p4] "=&v"(p4)
      : [x0] "v"(x.v[0]), [x1] "v"(x.v[1]), [x2] "v"(x.v[2]),
        [x3] "v"(x.v[3]), [x4] "v"(x.v[4]),
        [y0] "v"(y.v[0]), [y1] "v"(y.v[1]), [y2] "v"(y.v[2]),
        [y3] "v"(y.v[3]), [y4] "v"(y.v[4]),
        [mod0] "s"(p().v[0]), [mod1] "s"(p().v[1]), [mod2] "s"(p().v[2]),
        [mod3] "s"(p().v[3]), [mod4] "s"(p().v[4]));
  return z;
#else

  // Portable fallback for NVIDIA, older AMD GPUs, and wave64 compilation.
  // A 128-bit intermediate preserves borrow even when y's limb plus borrow wraps.
  uint128_t diff = (uint128_t)x.w.v0 - y.w.v0;
  z.w.v0 = (uint64_t)diff;
  diff = (uint128_t)x.w.v1 - y.w.v1 - ((diff >> 64) & 1);
  z.w.v1 = (uint64_t)diff;
  z.w.v2 = x.w.v2 - y.w.v2 - ((diff >> 64) & 1);
  if(z.w.v2 & 0x08) {
    z = add_raw(z, p());
  }
  return z;
#endif
}

__device__ uint131_t Curve<CURVE_ID_ECP131>::add(uint131_t x, uint131_t y)
{
  uint131_t z = add_raw(x, y);

  // Reduce mod P
  if(is_less_than(p(), z)) {
    z = sub_raw(z, p());
  }
  return z;
}

// Final Montgomery reduction: 0 <= x < 2P.
__device__ uint131_t mod_p(uint131_t x)
{
#if defined(__HIP_DEVICE_COMPILE__) && (defined(__GFX11__) || defined(__GFX12__))
  if(__builtin_amdgcn_wavefrontsize() == 32) {
    uint131_t z;
    uint32_t borrow;
    // Select x on underflow, otherwise x - P (including x == P).
    asm("v_sub_co_u32 %[z0], %[cc], %[x0], %[p0]\n\t"
        "v_sub_co_ci_u32 %[z1], %[cc], %[x1], %[p1], %[cc]\n\t"
        "v_sub_co_ci_u32 %[z2], %[cc], %[x2], %[p2], %[cc]\n\t"
        "v_sub_co_ci_u32 %[z3], %[cc], %[x3], %[p3], %[cc]\n\t"
        "v_sub_co_ci_u32 %[z4], %[cc], %[x4], %[p4], %[cc]\n\t"
        "v_cndmask_b32 %[z0], %[z0], %[x0], %[cc]\n\t"
        "v_cndmask_b32 %[z1], %[z1], %[x1], %[cc]\n\t"
        "v_cndmask_b32 %[z2], %[z2], %[x2], %[cc]\n\t"
        "v_cndmask_b32 %[z3], %[z3], %[x3], %[cc]\n\t"
        "v_cndmask_b32 %[z4], %[z4], %[x4], %[cc]"
        : [z0] "=&v"(z.v[0]), [z1] "=&v"(z.v[1]), [z2] "=&v"(z.v[2]),
          [z3] "=&v"(z.v[3]), [z4] "=&v"(z.v[4]), [cc] "=&s"(borrow)
        : [x0] "v"(x.v[0]), [x1] "v"(x.v[1]), [x2] "v"(x.v[2]),
          [x3] "v"(x.v[3]), [x4] "v"(x.v[4]),
          [p0] "s"(_p131_p.v[0]), [p1] "s"(_p131_p.v[1]), [p2] "s"(_p131_p.v[2]),
          [p3] "s"(_p131_p.v[3]), [p4] "s"(_p131_p.v[4]));
    return z;
  }
#endif
  if(!is_less_than(x, _p131_p)) {
    x = sub_raw(x, _p131_p);
  }
  return x;
}

// CIOS Montgomery multiplication for P131
// n=5 limbs of 32 bits, R=2^160, mp = -p^{-1} mod 2^32
// p.v[4] = 4 is hardcoded as a left-shift to save one multiply per iteration.
__device__ uint131_t Curve<CURVE_ID_ECP131>::mul(uint131_t a, uint131_t b)
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

  // REDC gives result < P + P*P/R. Most results are below P, so compare
  // the upper 35 bits first and only run the full subtraction near P.
  const uint64_t upper = ((uint64_t)result.v[4] << 32) | result.v[3];
  if(__builtin_expect(upper >= 0x48e1d43f2ULL, 0)) {
    result = mod_p(result);
  }
  return result;
}

// Squaring needs only 15 distinct limb products instead of multiplication's 25.
// The doubled cross terms need more than 64 bits before each carry is extracted.
// Follow the symmetric product with Montgomery REDC, R = 2^160.
__device__ uint131_t Curve<CURVE_ID_ECP131>::square(uint131_t a)
{
  uint64_t t[10];
  uint128_t c = 0;
  c += (uint128_t)((uint64_t)a.v[0] * a.v[0]);
  t[0] = (uint32_t)c;
  c >>= 32;
  c += ((uint128_t)((uint64_t)a.v[0] * a.v[1]) << 1);
  t[1] = (uint32_t)c;
  c >>= 32;
  c += ((uint128_t)((uint64_t)a.v[0] * a.v[2]) << 1) + (uint128_t)((uint64_t)a.v[1] * a.v[1]);
  t[2] = (uint32_t)c;
  c >>= 32;
  c += ((uint128_t)((uint64_t)a.v[0] * a.v[3]) << 1) + ((uint128_t)((uint64_t)a.v[1] * a.v[2]) << 1);
  t[3] = (uint32_t)c;
  c >>= 32;
  c += ((uint128_t)((uint64_t)a.v[0] * a.v[4]) << 1) +
       ((uint128_t)((uint64_t)a.v[1] * a.v[3]) << 1) +
       (uint128_t)((uint64_t)a.v[2] * a.v[2]);
  t[4] = (uint32_t)c;
  c >>= 32;
  c += ((uint128_t)((uint64_t)a.v[1] * a.v[4]) << 1) + ((uint128_t)((uint64_t)a.v[2] * a.v[3]) << 1);
  t[5] = (uint32_t)c;
  c >>= 32;
  c += ((uint128_t)((uint64_t)a.v[2] * a.v[4]) << 1) + (uint128_t)((uint64_t)a.v[3] * a.v[3]);
  t[6] = (uint32_t)c;
  c >>= 32;
  c += ((uint128_t)((uint64_t)a.v[3] * a.v[4]) << 1);
  t[7] = (uint32_t)c;
  c >>= 32;
  c += (uint128_t)((uint64_t)a.v[4] * a.v[4]);
  t[8] = (uint32_t)c;
  c >>= 32;
  t[9] = (uint32_t)c;

  const uint32_t mp = _p131_k.v[0];
  const uint32_t p0 = _p131_p.v[0];
  const uint32_t p1 = _p131_p.v[1];
  const uint32_t p2 = _p131_p.v[2];
  const uint32_t p3 = _p131_p.v[3];
#pragma unroll
  for(int i = 0; i < 5; ++i) {
    const uint32_t m = (uint32_t)t[i] * mp;
    uint64_t sum = (uint64_t)m * p0 + t[i];
    uint64_t carry = sum >> 32;
    sum = (uint64_t)m * p1 + t[i + 1] + carry;
    t[i + 1] = (uint32_t)sum;
    carry = sum >> 32;
    sum = (uint64_t)m * p2 + t[i + 2] + carry;
    t[i + 2] = (uint32_t)sum;
    carry = sum >> 32;
    sum = (uint64_t)m * p3 + t[i + 3] + carry;
    t[i + 3] = (uint32_t)sum;
    carry = sum >> 32;
    sum = ((uint64_t)m << 2) + t[i + 4] + carry;
    t[i + 4] = (uint32_t)sum;
    // Defer this carry: the next iteration normalizes it as t[i + 4].
    t[i + 5] += sum >> 32;
  }
  uint131_t result;
#pragma unroll
  for(int i = 0; i < 5; ++i) {
    result.v[i] = (uint32_t)t[i + 5];
  }
  return mod_p(result);
}

#endif
