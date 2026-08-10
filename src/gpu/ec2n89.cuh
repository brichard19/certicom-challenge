#ifndef _EC2N89_CUH
#define _EC2N89_CUH

#include "math_common.cuh"
#include "shared_types.h"

__constant__ uint131_t _ec2n89_a = {{0x660b75e77315e94e, 0x000000000095aa3e, 0x00000000}};
__constant__ uint131_t _ec2n89_b = {{0xc6c54021d1bf0a72, 0x0000000001ac2701, 0x00000000}};

template <> struct Curve<CURVE_ID_EC2N89> {
  __device__ static uint131_t one()
  {
    uint131_t value = {};
    value.v[0] = 1;
    return value;
  }
  __device__ static uint131_t a() { return _ec2n89_a; }
  __device__ static uint131_t b() { return _ec2n89_b; }

  __device__ static uint131_t add(uint131_t x, uint131_t y);
  __device__ static uint131_t sub(uint131_t x, uint131_t y);
  __device__ static uint131_t mul(uint131_t x, uint131_t y);
  __device__ static uint131_t square(uint131_t x);
  __device__ static uint131_t inv(uint131_t x);
};

__device__ uint131_t Curve<CURVE_ID_EC2N89>::add(uint131_t x, uint131_t y)
{
  uint131_t z;
  for(int i = 0; i < 5; i++) {
    z.v[i] = x.v[i] ^ y.v[i];
  }
  return z;
}

__device__ uint131_t Curve<CURVE_ID_EC2N89>::sub(uint131_t x, uint131_t y) { return add(x, y); }

__device__ uint131_t Curve<CURVE_ID_EC2N89>::mul(uint131_t x, uint131_t y)
{
  uint131_t product = {};

  for(int i = 0; i < 89; i++) {
    int bit = i % 32;
    int word = i / 32;

    if(y.v[word] & (uint32_t(1) << bit)) {
      product = add(product, x);
    }

    x.v[2] = (x.v[2] << 1) | (x.v[1] >> 31);
    x.v[1] = (x.v[1] << 1) | (x.v[0] >> 31);
    x.v[0] <<= 1;

    // x^89 = x^38 + 1 modulo x^89 + x^38 + 1.
    if(x.v[2] & 0x02000000) {
      x.v[0] ^= 0x00000001;
      x.v[1] ^= 0x00000040;
      x.v[2] ^= 0x02000000;
    }
  }

  return product;
}

__device__ uint64_t spread_bits_32_ec2n89(uint32_t x)
{
  uint64_t z = x;
  z = (z | (z << 16)) & 0x0000ffff0000ffffULL;
  z = (z | (z << 8)) & 0x00ff00ff00ff00ffULL;
  z = (z | (z << 4)) & 0x0f0f0f0f0f0f0f0fULL;
  z = (z | (z << 2)) & 0x3333333333333333ULL;
  z = (z | (z << 1)) & 0x5555555555555555ULL;
  return z;
}

__device__ uint131_t Curve<CURVE_ID_EC2N89>::square(uint131_t x)
{
  uint64_t w0 = spread_bits_32_ec2n89(x.v[0]);
  uint64_t w1 = spread_bits_32_ec2n89(x.v[1]);
  uint64_t w2 = spread_bits_32_ec2n89(x.v[2]);

  // Split at degree 89 and reduce x^89*h to (x^38 + 1)*h.
  uint64_t h0 = (w1 >> 25) | (w2 << 39);
  uint64_t h1 = w2 >> 25;

  uint64_t r0 = w0 ^ h0 ^ (h0 << 38);
  uint64_t r1 = (w1 & 0x1ffffffULL) ^ h1 ^ (h0 >> 26) ^ (h1 << 38);

  // The first fold can leave terms x^89 through x^125. One more fold
  // places all of them below degree 89.
  uint64_t high = r1 >> 25;
  r1 &= 0x1ffffffULL;
  r0 ^= high ^ (high << 38);
  r1 ^= high >> 26;

  uint131_t result = {};
  result.w.v0 = r0;
  result.w.v1 = r1;
  return result;
}

__device__ uint131_t Curve<CURVE_ID_EC2N89>::inv(uint131_t x)
{
  // Build x^(2^k-1), then square x^(2^88-1) to obtain x^(2^89-2).
  uint131_t x2 = mul(square(x), x);

  uint131_t x4 = x2;
  for(int i = 0; i < 2; i++)
    x4 = square(x4);
  x4 = mul(x4, x2);

  uint131_t x8 = x4;
  for(int i = 0; i < 4; i++)
    x8 = square(x8);
  x8 = mul(x8, x4);

  uint131_t x16 = x8;
  for(int i = 0; i < 8; i++)
    x16 = square(x16);
  x16 = mul(x16, x8);

  uint131_t x32 = x16;
  for(int i = 0; i < 16; i++)
    x32 = square(x32);
  x32 = mul(x32, x16);

  uint131_t x64 = x32;
  for(int i = 0; i < 32; i++)
    x64 = square(x64);
  x64 = mul(x64, x32);

  uint131_t x80 = x64;
  for(int i = 0; i < 16; i++)
    x80 = square(x80);
  x80 = mul(x80, x16);

  uint131_t x88 = x80;
  for(int i = 0; i < 8; i++)
    x88 = square(x88);
  x88 = mul(x88, x8);

  return square(x88);
}

#endif
