#ifndef _P131_CUH
#define _P131_CUH

#include "uint131.cuh"

__constant__ uint131_t _p131_p = {{0x194c43186b3abc0b, 0x8e1d43f293469e33, 0x4}};
__constant__ uint131_t _p131_k = {{0xe0587d72985b105d, 0xf1fd54b0309e1ab9, 0x7cfd70cf}};
__constant__ uint131_t _p131_r2 = {{0xf95d709f92600513, 0xf3d6fa1fb65ef639, 0x3}};
__constant__ uint131_t _p131_one = {{0x6e7743da32b6d0c7, 0x88c614d64c1a8f0b, 0x0}};
__constant__ uint131_t _p131_a = {{0xe7f7f250cee8709a, 0xacd15fe1a8ec1522, 0x0}};
__constant__ uint131_t _p131_b = {{0xc85087e5ab4eca9e, 0xde124657d7ba5851, 0x2}};

template<> struct Curve<131> {
  __device__ static uint131_t p() { return _p131_p;};
  __device__ static uint131_t k() { return _p131_k;};
  __device__ static uint131_t r2() { return _p131_r2;};
  __device__ static uint131_t one() { return _p131_one;};
  __device__ static uint131_t a() { return _p131_a;};
  __device__ static uint131_t b() { return _p131_b;};

  __device__ static uint131_t sub(uint131_t x, uint131_t y);
  __device__ static uint131_t add(uint131_t x, uint131_t y);
  
  __device__ static uint262_t mul(uint131_t x, uint131_t y);
  __device__ static uint262_t square(uint131_t x);

  __device__ static uint131_t mul_shift_160(uint160_t a, uint131_t b);

  __device__ static uint131_t high_bits(uint262_t x);
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

__device__ uint131_t Curve<131>::high_bits(uint262_t x)
{
  uint131_t hi;

  hi.w.v0 = (x.v[2] >> 32) | (x.v[3] << 32);
  hi.w.v1 = (x.v[3] >> 32) | (x.v[4] << 32);
  hi.w.v2 = 0;

  return hi;
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

// 131 x 131 -> 262 multiplication
__device__ uint262_t Curve<131>::mul(uint131_t a, uint131_t b)
{
  uint262_t tmp;
  uint64_t high = 0;

  // a0 * b0
  uint128_t t = (uint128_t)a.w.v0 * b.w.v0;
  tmp.v[0] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a.w.v0 * b.w.v1 + high;
  tmp.v[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b2
  t = (uint128_t)a.w.v0 * b.w.v2 + high;
  tmp.v[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  tmp.v[3] = high;

  // a1 * b0
  t = (uint128_t)a.w.v1 * b.w.v0 + tmp.v[1];
  tmp.v[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * b1 
  t = (uint128_t)a.w.v1 * b.w.v1 + tmp.v[2] + high;
  tmp.v[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * b2
  t = (uint128_t)a.w.v1 * b.w.v2 + tmp.v[3] + high;
  tmp.v[3] = (uint64_t)t;
  uint32_t high32 = (uint32_t)(t >> 64);

  tmp.v[4] = high32;

  // a2 * b0
  t = (uint128_t)a.w.v2* b.w.v0 + tmp.v[2];
  tmp.v[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * b1 
  t = (uint128_t)a.w.v2 * b.w.v1 + tmp.v[3] + high;
  tmp.v[3] = (uint64_t)t;
  high32 = (uint32_t)(t >> 64);

  // a2 * b2
  // The final word is only at most 6 bits, so no 128-bit mul needed
  uint32_t t32 = (uint32_t)a.w.v2 * b.w.v2 + tmp.v[4] + high32;
  tmp.v[4] = t32;

  return tmp;
}

// 131 x 131 -> 262 multiplication
__device__ uint262_t Curve<131>::square(uint131_t a)
{
  uint262_t tmp;
  uint64_t high = 0;

  // a0 * a0
  uint128_t t = (uint128_t)a.w.v0 * a.w.v0;
  tmp.v[0] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * a1 
  t = (uint128_t)a.w.v0 * a.w.v1 + high;
  tmp.v[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * a2
  t = (uint128_t)a.w.v0 * a.w.v2 + high;
  tmp.v[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  tmp.v[3] = high;

  // a1 * a0
  t = (uint128_t)a.w.v1 * a.w.v0 + tmp.v[1];
  tmp.v[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * a1 
  t = (uint128_t)a.w.v1 * a.w.v1 + tmp.v[2] + high;
  tmp.v[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * a2
  t = (uint128_t)a.w.v1 * a.w.v2 + tmp.v[3] + high;
  tmp.v[3] = (uint64_t)t;
  uint32_t high32 = (uint32_t)(t >> 64);

  tmp.v[4] = high32;

  // a2 * a0
  t = (uint128_t)a.w.v2 * a.w.v0 + tmp.v[2];
  tmp.v[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * a1 
  t = (uint128_t)a.w.v2 * a.w.v1 + tmp.v[3] + high;
  tmp.v[3] = (uint64_t)t;
  high32 = (uint32_t)(t >> 64);

  // a2 * a2
  // The final word is only at most 6 bits, so no 128-bit mul needed
  uint32_t t32 = (uint32_t)a.w.v2 * a.w.v2 + tmp.v[4] + high32;
  tmp.v[4] = t32;

  return tmp;
}

__device__ uint131_t Curve<131>::mul_shift_160(uint160_t a, uint131_t b)
{
  uint64_t tmp[5];
  uint64_t high = 0;

  // a0 * b0
  uint128_t t = (uint128_t)a.w.v0 * b.w.v0;
  high = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a.w.v0 * b.w.v1 + high;
  tmp[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b2
  t = (uint128_t)a.w.v0 * b.w.v2 + high;
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  tmp[3] = high;

  // a1 * b0
  t = (uint128_t)a.w.v1 * b.w.v0 + tmp[1];
  tmp[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);
  
  // a1 * b1 
  t = (uint128_t)a.w.v1 * b.w.v1 + tmp[2] + high;
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);
  
  // a1 * b2
  t = (uint128_t)a.w.v1 * b.w.v2 + tmp[3] + high;
  tmp[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);
  
  tmp[4] = high;
  
  // a2 * b0
  t = (uint128_t)a.w.v2 * b.w.v0 + tmp[2];
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);
  
  // a2 * b1 
  t = (uint128_t)a.w.v2 * b.w.v1 + tmp[3] + high;
  tmp[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);
  
  // a2 * b2
  // The final word is only at most 35 bits, so no 128-bit mul needed
  uint64_t t64 = (uint64_t)a.w.v2 * b.w.v2 + tmp[4] + high;
  tmp[4] = t64;

  uint131_t product;
  // Divide by 2^160
  product.w.v0 = (tmp[2] >> 32) | (tmp[3] << 32);
  product.w.v1 = (tmp[3] >> 32) | (tmp[4] << 32);
  product.w.v2 = (tmp[4] >> 32);

  return product;
}

#endif