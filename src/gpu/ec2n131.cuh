#ifndef _EC2N131_CUH
#define _EC2N131_CUH

#include "math_common.cuh"
#include "shared_types.h"

//__constant__ uint131_t _ec2n131_f = {{0x0000000000002007, 0x0000000000000000, 0x00000008}};
//__constant__ uint131_t _ec2n131_one = {{0x0000000000000001, 0x0000000000000000, 0x00000000}};
__constant__ uint131_t _ec2n131_a = {{0xa1a14f2c9e44352e, 0xebcb7eecc296a1c4, 0x00000007}};
__constant__ uint131_t _ec2n131_b = {{0x0093bdd622a61d81, 0x610b0a57c73649ad, 0x00000000}};

template <> struct Curve<CURVE_ID_EC2N131> {
  __device__ static uint131_t one()
  {
    uint131_t value = {};
    value.v[0] = 1;
    return value;
  }
  __device__ static uint131_t a() { return _ec2n131_a; }
  __device__ static uint131_t b() { return _ec2n131_b; }

  __device__ static uint131_t add(uint131_t x, uint131_t y);
  __device__ static uint131_t sub(uint131_t x, uint131_t y);
  __device__ static uint131_t mul(uint131_t x, uint131_t y);
  __device__ static uint131_t square(uint131_t x);
  __device__ static uint131_t inv(uint131_t x);
};

__device__ uint131_t Curve<CURVE_ID_EC2N131>::add(uint131_t x, uint131_t y)
{
  uint131_t z;
  for(int i = 0; i < 5; i++) {
    z.v[i] = x.v[i] ^ y.v[i];
  }
  return z;
}

__device__ uint131_t Curve<CURVE_ID_EC2N131>::sub(uint131_t x, uint131_t y) { return add(x, y); }

__device__ uint131_t Curve<CURVE_ID_EC2N131>::mul(uint131_t x, uint131_t y)
{
  uint131_t product = {};

  for(int i = 0; i < 131; i++) {
    int bit = i % 32;
    int word = i / 32;

    if(y.v[word] & (uint32_t(1) << bit)) {
      product = add(product, x);
    }

    x.v[4] = (x.v[4] << 1) | (x.v[3] >> 31);
    x.v[3] = (x.v[3] << 1) | (x.v[2] >> 31);
    x.v[2] = (x.v[2] << 1) | (x.v[1] >> 31);
    x.v[1] = (x.v[1] << 1) | (x.v[0] >> 31);
    x.v[0] <<= 1;

    if(x.v[4] & 0x08) {
      // x = add(x, _ec2n131_f);
      x.v[0] ^= 0x2007;
      x.v[4] ^= 0x08;
    }
  }

  return product;
}

__device__ uint64_t spread_bits_32(uint32_t x)
{
  uint64_t z = x;
  z = (z | (z << 16)) & 0x0000ffff0000ffffULL;
  z = (z | (z << 8)) & 0x00ff00ff00ff00ffULL;
  z = (z | (z << 4)) & 0x0f0f0f0f0f0f0f0fULL;
  z = (z | (z << 2)) & 0x3333333333333333ULL;
  z = (z | (z << 1)) & 0x5555555555555555ULL;
  return z;
}

__device__ uint131_t Curve<CURVE_ID_EC2N131>::square(uint131_t x)
{
  // Squaring inserts a zero between every pair of coefficients. Each input
  // limb therefore expands independently into one 64-bit word.
  uint64_t w0 = spread_bits_32(x.v[0]);
  uint64_t w1 = spread_bits_32(x.v[1]);
  uint64_t w2 = spread_bits_32(x.v[2]);
  uint64_t w3 = spread_bits_32(x.v[3]);
  uint64_t w4 = spread_bits_32(x.v[4]);

  // Split the 261-bit square at degree 131. For high polynomial h,
  // x^131*h reduces to (x^13 + x^2 + x + 1)*h.
  uint64_t h0 = (w2 >> 3) | (w3 << 61);
  uint64_t h1 = (w3 >> 3) | (w4 << 61);
  uint64_t h2 = w4 >> 3;

  uint64_t r0 = w0 ^ h0 ^ (h0 << 1) ^ (h0 << 2) ^ (h0 << 13);
  uint64_t r1 = w1 ^ h1 ^ (h1 << 1) ^ (h0 >> 63) ^ (h1 << 2) ^ (h0 >> 62) ^ (h1 << 13) ^ (h0 >> 51);
  uint64_t r2 =
      (w2 & 0x7) ^ h2 ^ (h2 << 1) ^ (h1 >> 63) ^ (h2 << 2) ^ (h1 >> 62) ^ (h2 << 13) ^ (h1 >> 51);

  // The first fold can only leave terms x^131 through x^142. Fold those
  // once more; their reductions all fit in the low 64-bit word.
  uint64_t high = r2 >> 3;
  r2 &= 0x7;
  r0 ^= high ^ (high << 1) ^ (high << 2) ^ (high << 13);

  uint131_t result;
  result.w.v0 = r0;
  result.w.v1 = r1;
  result.w.v2 = (uint32_t)r2;
  return result;
}

__device__ uint131_t Curve<CURVE_ID_EC2N131>::inv(uint131_t x)
{
  // Addition chain for 2^131 - 2. Exponent shifts are repeated squaring;
  // exponent additions are field multiplication.
  uint131_t x10 = square(x);
  uint131_t x11 = mul(x, x10);

  uint131_t x1100 = x11;
  for(int i = 0; i < 2; i++) {
    x1100 = square(x1100);
  }

  uint131_t x1111 = mul(x11, x1100);

  uint131_t x11110000 = x1111;
  for(int i = 0; i < 4; i++) {
    x11110000 = square(x11110000);
  }

  uint131_t x11111111 = mul(x1111, x11110000);

  uint131_t x16 = x11111111;
  for(int i = 0; i < 8; i++) {
    x16 = square(x16);
  }
  x16 = mul(x16, x11111111);

  uint131_t x32 = x16;
  for(int i = 0; i < 16; i++) {
    x32 = square(x32);
  }
  x32 = mul(x32, x16);

  uint131_t x64 = x32;
  for(int i = 0; i < 32; i++) {
    x64 = square(x64);
  }
  x64 = mul(x64, x32);

  uint131_t x128 = x64;
  for(int i = 0; i < 64; i++) {
    x128 = square(x128);
  }
  x128 = mul(x128, x64);

  uint131_t x130 = square(x128);
  x130 = square(x130);
  x130 = mul(x130, x11);

  return square(x130);
}

#endif
