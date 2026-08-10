#ifndef _EC2N79_CUH
#define _EC2N79_CUH

#include "math_common.cuh"
#include "shared_types.h"

__constant__ uint131_t _ec2n79_a = {{0x38a8f66d7f4c385f, 0x0000000000004a2e, 0x00000000}};
__constant__ uint131_t _ec2n79_b = {{0xb31c6becc03d68a7, 0x0000000000002c0b, 0x00000000}};

template <> struct Curve<CURVE_ID_EC2N79> {
  __device__ static uint131_t one()
  {
    uint131_t value = {};
    value.v[0] = 1;
    return value;
  }
  __device__ static uint131_t a() { return _ec2n79_a; }
  __device__ static uint131_t b() { return _ec2n79_b; }

  __device__ static uint131_t add(uint131_t x, uint131_t y);
  __device__ static uint131_t sub(uint131_t x, uint131_t y);
  __device__ static uint131_t mul(uint131_t x, uint131_t y);
  __device__ static uint131_t square(uint131_t x);
  __device__ static uint131_t inv(uint131_t x);
};

__device__ uint131_t Curve<CURVE_ID_EC2N79>::add(uint131_t x, uint131_t y)
{
  uint131_t z;
  for(int i = 0; i < 5; i++) {
    z.v[i] = x.v[i] ^ y.v[i];
  }
  return z;
}

__device__ uint131_t Curve<CURVE_ID_EC2N79>::sub(uint131_t x, uint131_t y) { return add(x, y); }

__device__ uint131_t Curve<CURVE_ID_EC2N79>::mul(uint131_t x, uint131_t y)
{
  uint131_t product = {};

  for(int i = 0; i < 79; i++) {
    int bit = i % 32;
    int word = i / 32;

    if(y.v[word] & (uint32_t(1) << bit)) {
      product = add(product, x);
    }

    x.v[2] = (x.v[2] << 1) | (x.v[1] >> 31);
    x.v[1] = (x.v[1] << 1) | (x.v[0] >> 31);
    x.v[0] <<= 1;

    // x^79 = x^9 + 1 modulo x^79 + x^9 + 1.
    if(x.v[2] & 0x00008000) {
      x.v[0] ^= 0x00000201;
      x.v[2] ^= 0x00008000;
    }
  }

  return product;
}

__device__ uint64_t spread_bits_32_ec2n79(uint32_t x)
{
  uint64_t z = x;
  z = (z | (z << 16)) & 0x0000ffff0000ffffULL;
  z = (z | (z << 8)) & 0x00ff00ff00ff00ffULL;
  z = (z | (z << 4)) & 0x0f0f0f0f0f0f0f0fULL;
  z = (z | (z << 2)) & 0x3333333333333333ULL;
  z = (z | (z << 1)) & 0x5555555555555555ULL;
  return z;
}

__device__ uint131_t Curve<CURVE_ID_EC2N79>::square(uint131_t x)
{
  uint64_t w0 = spread_bits_32_ec2n79(x.v[0]);
  uint64_t w1 = spread_bits_32_ec2n79(x.v[1]);
  uint64_t w2 = spread_bits_32_ec2n79(x.v[2]);

  // Split at degree 79 and reduce x^79*h to (x^9 + 1)*h.
  uint64_t h0 = (w1 >> 15) | (w2 << 49);
  uint64_t h1 = w2 >> 15;

  uint64_t r0 = w0 ^ h0 ^ (h0 << 9);
  uint64_t r1 = (w1 & 0x7fffULL) ^ h1 ^ (h0 >> 55) ^ (h1 << 9);

  // The first fold can leave terms x^79 through x^86. Fold them once more.
  uint64_t high = r1 >> 15;
  r1 &= 0x7fffULL;
  r0 ^= high ^ (high << 9);

  uint131_t result = {};
  result.w.v0 = r0;
  result.w.v1 = r1;
  return result;
}

__device__ uint131_t Curve<CURVE_ID_EC2N79>::inv(uint131_t x)
{
  // Build x^(2^k-1), then square x^(2^78-1) to obtain x^(2^79-2).
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

  uint131_t x72 = x64;
  for(int i = 0; i < 8; i++)
    x72 = square(x72);
  x72 = mul(x72, x8);

  uint131_t x76 = x72;
  for(int i = 0; i < 4; i++)
    x76 = square(x76);
  x76 = mul(x76, x4);

  uint131_t x78 = x76;
  for(int i = 0; i < 2; i++)
    x78 = square(x78);
  x78 = mul(x78, x2);

  return square(x78);
}

#endif
