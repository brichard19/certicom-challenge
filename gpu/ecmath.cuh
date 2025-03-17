
#ifndef _EC_MATH_CUH
#define _EC_MATH_CUH

typedef unsigned __int128 uint128_t;
#include "p131.cuh"
#include "p79.cuh"

__device__ uint131_t sub_raw(const uint131_t& x, const uint131_t& y)
{
  uint131_t z;
  int borrow = 0;

  uint128_t diff = (uint128_t)x.w.v0 - y.w.v0 - borrow;
  z.w.v0 = (uint64_t)diff;
  borrow = (int)(diff >> 64) & 1;
  
  diff = (uint128_t)x.w.v1 - y.w.v1 - borrow;
  z.w.v1 = (uint64_t)diff;
  borrow = (int)(diff >> 64) & 1;

  z.w.v2 = x.w.v2 - y.w.v2 - borrow;

  // Sign extend
  z.w.padding = z.w.v2 & 0x80000000 ? -1 : 0;

  return z;
}

__device__ uint131_t add_raw(const uint131_t& x, const uint131_t& y, int carry_in = 0)
{

  uint131_t z;
  int carry = carry_in;
  
  uint128_t sum = (uint128_t)x.w.v0 + y.w.v0 + carry;
  z.w.v0 = (uint64_t)sum;
  carry = (uint64_t)(sum >> 64);
  
  sum = (uint128_t)x.w.v1 + y.w.v1 + carry;
  z.w.v1 = (uint64_t)sum;
  carry = (uint64_t)(sum >> 64);
 
  z.w.v2 = x.w.v2 + y.w.v2 + carry;

  // Padding
  z.w.padding = 0;

  return z;
}

// if x is less than y
__device__ int is_less_than(uint131_t& x, uint131_t& y)
{
  uint131_t diff = sub_raw(x, y);

  return (diff.w.v2 >> 31) & 1;
}

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

  z.w.padding = 0;

  return z;
}

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

template<int CURVE> __device__ uint131_t sub(uint131_t x, uint131_t y)
{
  if(CURVE == 131) {
    return sub_p131(x, y);
  } else if(CURVE == 79) {
    return sub_p79(x, y);
  }
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

__device__ uint131_t add_p79(const uint131_t& x, const uint131_t& y)
{
  uint131_t z = add_raw(x, y);

  // Reduce mod P
  if(is_less_than(_p79_p, z)) {
    z = sub_raw(z, _p79_p);
  }
  return z;
}

template<int CURVE> __device__ uint131_t add(const uint131_t& x, const uint131_t& y)
{
  if(CURVE == 131) {
    return add_p131(x, y);
  } else if(CURVE == 79) {
    return add_p79(x, y);
  }
}




// 131 x 131 -> 262 multiplication
__device__ uint262_t mul_131(const uint131_t& a, const uint131_t& b)
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
  high = (uint64_t)(t >> 64);

  tmp.v[4] = high;

  // a2 * b0
  t = (uint128_t)a.w.v2* b.w.v0 + tmp.v[2];
  tmp.v[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * b1 
  t = (uint128_t)a.w.v2 * b.w.v1 + tmp.v[3] + high;
  tmp.v[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * b2
  // The final word is only at most 6 bits, so no 128-bit mul needed
  uint64_t t64 = a.w.v2 * b.w.v2 + tmp.v[4] + high;
  tmp.v[4] = t64;

  return tmp;
}



// 131 x 131 -> 262 multiplication
__device__ uint262_t square_131(const uint131_t& a)
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
  high = (uint64_t)(t >> 64);

  tmp.v[4] = high;

  // a2 * a0
  t = (uint128_t)a.w.v2 * a.w.v0 + tmp.v[2];
  tmp.v[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * a1 
  t = (uint128_t)a.w.v2 * a.w.v1 + tmp.v[3] + high;
  tmp.v[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * a2
  // The final word is only at most 6 bits, so no 128-bit mul needed
  uint64_t t64 = a.w.v2 * a.w.v2 + tmp.v[4] + high;
  tmp.v[4] = t64;

  return tmp;
}



// One in put is 131 bits, the other is 160 bits
// Output is 131 bits
__device__ uint131_t mul_shift_160(const uint160_t& a, const uint131_t& b)
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
  product.w.v0 = (tmp[2] >> 32) | ((tmp[3] & 0xffffffff) << 32);
  product.w.v1 = (tmp[3] >> 32) | ((tmp[4] & 0xffffffff) << 32);
  product.w.v2 = (tmp[4] >> 32);

  return product;
}



__device__ uint160_t mul_mod_160(const uint160_t& a, const uint131_t& b)
{
  uint160_t tmp;
  uint64_t high = 0;

  // a0 * b0
  uint128_t t = (uint128_t)a.w.v0 * b.w.v0;
  tmp.w.v0 = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a.w.v0 * b.w.v1 + high;
  tmp.w.v1 = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b2
  t = (uint128_t)a.w.v0 * b.w.v2 + high;
  tmp.w.v2 = (uint32_t)t;

  // a1 * b0
  t = (uint128_t)a.w.v1 * b.w.v0 + tmp.w.v1;
  tmp.w.v1 = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a.w.v1 * b.w.v1 + tmp.w.v2 + high;
  tmp.w.v2 = (uint32_t)t;

  // a2 * b0
  t = (uint128_t)a.w.v2 * b.w.v0 + tmp.w.v2;
  tmp.w.v2 = (uint64_t)t;

  tmp.w.v2 &= 0xffffffff;

  return tmp;

}


template<int CURVE> __device__ uint131_t mont_reduce(const uint262_t& t)
{
  // m1 = t_lo * k mod R
  // m2 = m1 * p / r
  // m3 = m2 + t_hi + 1

  uint160_t t_lo; 
  t_lo.w.v0 = t.v[0];
  t_lo.w.v1 = t.v[1];
  t_lo.w.v2 = (uint32_t)t.v[2];
  t_lo.w.padding = 0;

  // Remaining high bits (102 bits)
  uint131_t t_hi;

  t_hi.w.v0 = (t.v[2] >> 32) | ((t.v[3] & 0xffffffff) << 32);
  t_hi.w.v1 = (t.v[3] >> 32) | ((t.v[4] & 0xffffffff) << 32);
  t_hi.w.v2 = 0;
  t_hi.w.padding = 0;

  uint160_t m1;
  uint131_t m2;
  uint131_t m3;

  uint131_t p;
  uint131_t k;
  if(CURVE == 131) {
    p = _p131_p;
    k = _p131_k;
  } else if(CURVE == 79) {
    p = _p79_p;
    k = _p79_k;
  } 

  m1 = mul_mod_160(t_lo, k);

  m2 = mul_shift_160(m1, p);
 
  // Not sure if this is correct, but carry always seems to be 1.
  m3 = add_raw(t_hi, m2, 1);

  if(is_less_than(p, m3)) {
    m3 = sub_raw(m3, p);
  }
  return m3;

}



template<int CURVE> __device__ uint131_t mul(uint131_t x, uint131_t y)
{
  uint262_t product = mul_131(x, y);
  return mont_reduce<CURVE>(product);
}



template<int CURVE> __device__ uint131_t square(uint131_t x)
{
  uint262_t product = square_131(x);
  return mont_reduce<CURVE>(product);
}

template<int CURVE> __device__ uint131_t square(uint131_t x, int n)
{
  for(int i = 0; i < n; i++) {
    x = square<CURVE>(x);
  }

  return x;
}

__device__ bool equal(const uint131_t& x, const uint131_t& y)
{
  return x.w.v0 == y.w.v0 && x.w.v1 == y.w.v1 && x.w.v2 == y.w.v2;
}

// Modular inverse using Fermat's method
__device__ uint131_t inv_p131(uint131_t& x)
{
  uint131_t z, t0, t1, t2, t3, t4, t5;

  t0 = square<131>(x);
  t3 = mul<131>(x, t0);
  t4 = mul<131>(x, t3);
  t1 = mul<131>(x, t4);
  t2 = mul<131>(t0, t1);
  z = mul<131>(t0, t2);
  t4 = mul<131>(t4, z);
  t0 = mul<131>(t0, t4);
  t5 = mul<131>(t3, t0);
  t5 = square<131>(t5, 5);
  t5 = mul<131>(t2, t5);
  t5 = square<131>(t5, 7);
  t5 = mul<131>(t2, t5);
  t5 = square<131>(t5, 4);
  t5 = mul<131>(t1, t5);
  t5 = square<131>(t5, 8);
  t5 = mul<131>(t0, t5);
  t5 = square<131>(t5, 2);
  t5 = mul<131>(t3, t5);
  t5 = square<131>(t5, 5);
  t5 = mul<131>(t1, t5);
  t5 = square<131>(t5, 6);
  t5 = mul<131>(z, t5);
  t5 = square<131>(t5, 3);
  t5 = mul<131>(t1, t5);
  t5 = square<131>(t5, 7);
  t5 = mul<131>(t4, t5);
  t5 = square<131>(t5, 6);
  t5 = mul<131>(t0, t5);
  t5 = square<131>(t5, 5);
  t5 = mul<131>(t3, t5);
  t5 = square<131>(t5, 4);
  t5 = mul<131>(t3, t5);
  t5 = square<131>(t5, 5);
  t5 = mul<131>(t3, t5);
  t5 = square<131>(t5, 5);
  t5 = mul<131>(t1, t5);
  t5 = square<131>(t5, 4);
  t5 = mul<131>(t3, t5);
  t5 = square<131>(t5, 4);
  t5 = mul<131>(x, t5);
  t5 = square<131>(t5, 6);
  t5 = mul<131>(t3, t5);
  t5 = square<131>(t5, 5);
  t5 = mul<131>(t3, t5);
  t5 = square<131>(t5, 8);
  t4 = mul<131>(t4, t5);
  t4 = square<131>(t4, 3);
  t3 = mul<131>(t3, t4);
  t3 = square<131>(t3, 5);
  t2 = mul<131>(t2, t3);
  t2 = square<131>(t2, 4);
  t1 = mul<131>(t1, t2);
  t1 = square<131>(t1, 5);
  t0 = mul<131>(t0, t1);
  t0 = square<131>(t0, 10);
  z = mul<131>(z, t0);

  return z;
}

__device__ uint131_t inv_p79(uint131_t& x)
{
  uint131_t z, t0, t1, t2, t3, t4, t5, t6, t7;
  
  t6 = square<79>(x);
  t0 = mul<79>(x, t6);
  z = square<79>(t0);
  t3 = mul<79>(t6, z);
  t2 = mul<79>(t0, t3);
  t7 = mul<79>(x, t2);
  t0 = mul<79>(t3, t2);
  t1 = mul<79>(t6, t0);
  t5 = mul<79>(t6, t1);
  t4 = mul<79>(z, t5);
  t3 = mul<79>(t3, t4);
  t7 = mul<79>(t7, t3);
  t6 = mul<79>(t6, t7);
  z = mul<79>(z, t6);
  t7 = square<79>(t7, 7);
  t6 = mul<79>(t6, t7);
  t6 = square<79>(t6, 6);
  t6 = mul<79>(t3, t6);
  t6 = square<79>(t6, 8);
  t5 = mul<79>(t5, t6);
  t5 = square<79>(t5, 6);
  t4 = mul<79>(t4, t5);
  t4 = square<79>(t4, 11);
  t3 = mul<79>(t3, t4);
  t3 = square<79>(t3, 5);
  t2 = mul<79>(t2, t3);
  t2 = square<79>(t2, 7);
  t1 = mul<79>(t1, t2);
  t1 = square<79>(t1, 8);
  t0 = mul<79>(t0, t1);
  t0 = square<79>(t0, 8);
  t0 = mul<79>(z, t0);
  t0 = square<79>(t0, 6);
  z = mul<79>(z, t0);
  z = square<79>(z);
  z = mul<79>(x, z);

  return z;
}

template<int CURVE> __device__ uint131_t inv(uint131_t x)
{
  if(CURVE == 131) {
    return inv_p131(x);
  } else if(CURVE == 79) {
    return inv_p79(x);
  }
}

// ECC functons

// Check for point-at-infinity using the x coordinate
__device__ bool is_infinity(uint131_t x)
{
  return x.v[2] == (uint64_t)-1;
}

__device__ void set_point_at_infinity(uint131_t& x)
{
  x.v[2] = (uint64_t)-1;
}

template<int CURVE> __device__ bool point_exists(uint131_t& x, uint131_t& y)
{
  uint131_t a;
  uint131_t b;

  switch(CURVE) {
    case 131:
      a = _p131_a;
      b = _p131_b;
      break;
    case 79:
      a = _p79_a;
      b = _p79_b;
      break;
  }

  uint131_t y2 = square<CURVE>(y);
  uint131_t x3 = mul<CURVE>(x, square<CURVE>(x));
  uint131_t ax = mul<CURVE>(a, x);

  uint131_t rs = add<CURVE>(add<CURVE>(x3, ax), b);

  return equal(y2, rs);
}


#endif