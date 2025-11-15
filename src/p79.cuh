#ifndef _P79_CUH
#define _P79_CUH

#include "uint131.cuh"

__constant__ uint131_t _p79_p = {{0x5177412aca899cf5, 0x62ce, 0x0}};
__constant__ uint131_t _p79_k = {{0x6e3655426732d0a3, 0xcafea9fd045a89b6, 0xe6bb05ec}};
__constant__ uint131_t _p79_r2 = {{0x7b0baef57de52417, 0xe79, 0x0}};
__constant__ uint131_t _p79_one = {{0x5447aa703f6abc5f, 0x1358, 0x0}};
__constant__ uint131_t _p79_a = {{0x732c9b460e3c3d, 0x1bb7, 0x0}};
__constant__ uint131_t _p79_b = {{0xc88edfd7d5b44610, 0x250c, 0x0}};

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


#endif