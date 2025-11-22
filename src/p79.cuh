#ifndef _P79_CUH
#define _P79_CUH

#include "uint131.cuh"

__constant__ uint131_t _p79_p = {{0x5177412aca899cf5, 0x62ce, 0x0}};
__constant__ uint131_t _p79_k = {{0x6e3655426732d0a3, 0xcafea9fd045a89b6, 0xe6bb05ec}};
__constant__ uint131_t _p79_r2 = {{0x7b0baef57de52417, 0xe79, 0x0}};
__constant__ uint131_t _p79_one = {{0x5447aa703f6abc5f, 0x1358, 0x0}};
__constant__ uint131_t _p79_a = {{0x732c9b460e3c3d, 0x1bb7, 0x0}};
__constant__ uint131_t _p79_b = {{0xc88edfd7d5b44610, 0x250c, 0x0}};


template<> struct Curve<79> {
  __device__ static uint131_t p() { return _p79_p;};
  __device__ static uint131_t k() { return _p79_k;};
  __device__ static uint131_t r2() { return _p79_r2;};
  __device__ static uint131_t one() { return _p79_one;};
  __device__ static uint131_t a() { return _p79_a;};
  __device__ static uint131_t b() { return _p79_b;};
};

__device__ uint131_t sub_p79(uint131_t x, uint131_t y)
{
  uint131_t z = {{0}};

  int borrow = 0;
  
  uint64_t diff = x.w.v0 - y.w.v0;
  z.w.v0 = diff;
  borrow = diff > x.w.v0 ? 1 : 0;
 
  // High word is only 15 bits
  diff = x.w.v1 - y.w.v1 - borrow;
  z.w.v1 = diff;
  borrow = diff & 0x8000;
 
  // Went below zero. Need to add P.
  if(borrow) {
    int carry = 0; 
    uint64_t sum = z.w.v0 + _p79_p.w.v0;
    z.w.v0 = sum;
    carry = sum < _p79_p.w.v0 ? 1 : 0;
    
    sum = z.w.v1 + _p79_p.w.v1 + carry;
    z.w.v1 = sum;
  }

  return z;
}

__device__ uint131_t add_p79(const uint131_t& x, const uint131_t& y)
{
  uint131_t z = add_raw(x, y);

  // Reduce mod P
  if(is_less_than(_p79_p, z)) {
    z = sub_raw(z, _p79_p);
  }
  return z;
}

// 131 x 131 -> 262 multiplication
__device__ uint262_t mul_p79(const uint131_t& a, const uint131_t& b)
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

  tmp.v[2] = high;

  tmp.v[3] = 0;

  // a1 * b0
  t = (uint128_t)a.w.v1 * b.w.v0 + tmp.v[1];
  tmp.v[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * b1 
  t = (uint128_t)a.w.v1 * b.w.v1 + tmp.v[2] + high;
  tmp.v[2] = (uint64_t)t;

  tmp.v[4] = 0;


  return tmp;
}


__device__ uint262_t square_p79(const uint131_t& a)
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

  tmp.v[2] = high;
  tmp.v[3] = 0;

  // a1 * a0
  t = (uint128_t)a.w.v1 * a.w.v0 + tmp.v[1];
  tmp.v[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * a1 
  t = (uint128_t)a.w.v1 * a.w.v1 + tmp.v[2] + high;
  tmp.v[2] = (uint64_t)t;

  tmp.v[4] = 0;

  return tmp;
}

__device__ uint131_t mul_shift_160_p79(const uint160_t& a, const uint131_t& b)
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

  tmp[2] = high;
  tmp[3] = 0;

  // a1 * b0
  t = (uint128_t)a.w.v1 * b.w.v0 + tmp[1];
  tmp[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);
  
  // a1 * b1 
  t = (uint128_t)a.w.v1 * b.w.v1 + tmp[2] + high;
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);
  tmp[3] = high;  

  // a2 * b0
  t = (uint128_t)a.w.v2 * b.w.v0 + tmp[2];
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);
  
  // a2 * b1 
  t = (uint128_t)a.w.v2 * b.w.v1 + tmp[3] + high;
  tmp[3] = (uint64_t)t;

  tmp[4] = 0;

  uint131_t product;
  // Divide by 2^160
  product.w.v0 = (tmp[2] >> 32) | (tmp[3] << 32);
  product.w.v1 = (tmp[3] >> 32) | (tmp[4] << 32);
  product.w.v2 = (tmp[4] >> 32);

  return product;
}

#endif