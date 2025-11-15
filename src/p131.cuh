#ifndef _P131_CUH
#define _P131_CUH

#include "uint131.cuh"

__constant__ uint131_t _p131_p = {{0x194c43186b3abc0b, 0x8e1d43f293469e33, 0x4}};
__constant__ uint131_t _p131_k = {{0xe0587d72985b105d, 0xf1fd54b0309e1ab9, 0x7cfd70cf}};
__constant__ uint131_t _p131_r2 = {{0xf95d709f92600513, 0xf3d6fa1fb65ef639, 0x3}};
__constant__ uint131_t _p131_one = {{0x6e7743da32b6d0c7, 0x88c614d64c1a8f0b, 0x0}};
__constant__ uint131_t _p131_a = {{0xe7f7f250cee8709a, 0xacd15fe1a8ec1522, 0x0}};
__constant__ uint131_t _p131_b = {{0xc85087e5ab4eca9e, 0xde124657d7ba5851, 0x2}};



__device__ uint131_t sub_p131(uint131_t x, uint131_t y)
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
    uint64_t sum = z.w.v0 + _p131_p.w.v0;
    z.w.v0 = sum;
    carry = sum < _p131_p.w.v0 ? 1 : 0;
    
    sum = z.w.v1 + _p131_p.w.v1 + carry;
    z.w.v1 = sum;
    carry = sum < _p131_p.w.v1 ? 1 : 0;
    
    z.w.v2 = z.w.v2 + _p131_p.w.v2 + carry;
  }


  return z;
}

__device__ uint131_t add_p131(const uint131_t& x, const uint131_t& y)
{
  uint131_t z = add_raw(x, y);

  // Reduce mod P
  if(is_less_than(_p131_p, z)) {
    z = sub_raw(z, _p131_p);
  }
  return z;
}


#endif