
#ifndef _EC_MATH_CUH
#define _EC_MATH_CUH

typedef unsigned __int128 uint128_t;

#if defined(CURVE_P131)
#include "p131.cuh"
#elif defined(CURVE_P79)
#include "p79.cuh"
#else
#error "Curve not defined"
#endif

__device__ uint131_t sub_p131(uint131_t x, uint131_t y)
{
  uint131_t z;

  int borrow = 0;
  
  uint64_t diff = x.v[0] - y.v[0];
  z.v[0] = diff;
  borrow = diff > x.v[0] ? 1 : 0;
  
  diff = x.v[1] - y.v[1] - borrow;
  z.v[1] = diff;
  borrow = diff > x.v[1] ? 1 : 0;
 
  // High word is only 3 bits
  diff = x.v[2] - y.v[2] - borrow;
  z.v[2] = diff;
  borrow = diff & 0x08;

  // Went below zero. Need to add P.
  if(borrow) {
    int carry = 0; 
    uint64_t sum = z.v[0] + _p.v[0];
    z.v[0] = sum;
    carry = sum < _p.v[0] ? 1 : 0;
    
    sum = z.v[1] + _p.v[1] + carry;
    z.v[1] = sum;
    carry = sum < _p.v[1] ? 1 : 0;
    
    z.v[2] = z.v[2] + _p.v[2] + carry;
  }

  return z;
}

__device__ uint131_t sub_p79(uint131_t x, uint131_t y)
{
  uint131_t z = {{0}};

  int borrow = 0;
  
  uint64_t diff = x.v[0] - y.v[0];
  z.v[0] = diff;
  borrow = diff > x.v[0] ? 1 : 0;
 
  // High word is only 15 bits
  diff = x.v[1] - y.v[1] - borrow;
  z.v[1] = diff;
  borrow = diff & 0x8000;
 
  // Went below zero. Need to add P.
  if(borrow) {
    int carry = 0; 
    uint64_t sum = z.v[0] + _p.v[0];
    z.v[0] = sum;
    carry = sum < _p.v[0] ? 1 : 0;
    
    sum = z.v[1] + _p.v[1] + carry;
    z.v[1] = sum;
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


__device__ uint131_t add(const uint131_t& x, const uint131_t& y)
{
  uint131_t z = {{0}};

  uint64_t carry = 0;

  for(int i = 0; i < 3; i++) {
    uint128_t sum = (uint128_t)x.v[i] + y.v[i] + carry;

    z.v[i] = (uint64_t)sum;
    carry = (uint64_t)(sum >> 64);
  }

  // Reduce mod P
  bool gte = true;
  for(int i = 2; i >= 0; i--) {
    if(z.v[i] > _p.v[i]) {
      break;
    } else if(z.v[i] < _p.v[i]) {
      gte = false;
      break;
    }
  }

  // Subtract P
  if(gte) {
    uint64_t borrow = 0;
    for(int i = 0; i < 3; i++) {
      uint128_t diff = (uint128_t)z.v[i] - _p.v[i] - borrow;
      borrow = (uint64_t)(diff >> 64) & 1;
      z.v[i] = (uint64_t)diff;
    }
  }

  return z;
}


// 160 x 160 -> 320 multiplication
__device__ void mul160x(const uint64_t* a, const uint64_t* b, uint64_t* product)
{
  uint64_t high = 0;

  // a0 * b0
  uint128_t t = (uint128_t)a[0] * b[0] + 0;
  product[0] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a[0] * b[1] + 0 + high;
  product[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b2
  t = (uint128_t)a[0] * b[2] + 0 + high;
  product[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  product[3] = high;

  // a1 * b0
  t = (uint128_t)a[1] * b[0] + product[1];
  product[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * b1 
  t = (uint128_t)a[1] * b[1] + product[2] + high;
  product[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * b2
  t = (uint128_t)a[1] * b[2] + product[3] + high;
  product[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  product[4] = high;

  // a2 * b0
  t = (uint128_t)a[2] * b[0] + product[2];
  product[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * b1 
  t = (uint128_t)a[2] * b[1] + product[3] + high;
  product[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * b2
  // The final word is only at most 6 bits, so no 128-bit mul needed
  uint64_t t64 = a[2] * b[2] + product[4] + high;
  product[4] = t64;
}

// 160 x 160 -> 320 squaring
__device__ void square160(const uint64_t* a, uint64_t* product)
{
  //uint64_t tmp[5] = {0};
  uint64_t high = 0;

  // a0 * a0
  uint128_t t = (uint128_t)a[0] * a[0] + 0;
  product[0] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * a1 
  t = (uint128_t)a[0] * a[1] + 0 + high;
  product[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * a2
  t = (uint128_t)a[0] * a[2] + 0 + high;
  product[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  product[3] = high;

  // a1 * a0
  t = (uint128_t)a[1] * a[0] + product[1];
  product[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * a1 
  t = (uint128_t)a[1] * a[1] + product[2] + high;
  product[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * a2
  t = (uint128_t)a[1] * a[2] + product[3] + high;
  product[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  product[4] = high;

  // a2 * a0
  t = (uint128_t)a[2] * a[0] + product[2];
  product[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * a1 
  t = (uint128_t)a[2] * a[1] + product[3] + high;
  product[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * a2
  // The final word is only at most 6 bits, so no 128-bit mul needed
  uint64_t t64 = a[2] * a[2] + product[4] + high;
  product[4] = t64;

}

__device__ void mont_mulDivR(const uint64_t* a, const uint64_t* b, uint64_t* product)
{
  uint64_t tmp[5] = {0};
  uint64_t high = 0;

  // a0 * b0
  uint128_t t = (uint128_t)a[0] * b[0] + tmp[0];
  //tmp[0] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a[0] * b[1] + tmp[1] + high;
  tmp[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b2
  t = (uint128_t)a[0] * b[2] + tmp[2] + high;
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  tmp[3] = high;

  // a1 * b0
  t = (uint128_t)a[1] * b[0] + tmp[1];
  tmp[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * b1 
  t = (uint128_t)a[1] * b[1] + tmp[2] + high;
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * b2
  t = (uint128_t)a[1] * b[2] + tmp[3] + high;
  tmp[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  tmp[4] = high;

  // a2 * b0
  t = (uint128_t)a[2] * b[0] + tmp[2];
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * b1 
  t = (uint128_t)a[2] * b[1] + tmp[3] + high;
  tmp[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * b2
  // The final word is only at most 6 bits, so no 128-bit mul needed
  uint64_t t64 = a[2] * b[2] + tmp[4] + high;
  tmp[4] = t64;

  product[0] = (tmp[2] >> 32) | ((tmp[3] & 0xffffffff) << 32);
  product[1] = (tmp[3] >> 32) | ((tmp[4] & 0xffffffff) << 32);
  product[2] = (tmp[4] >> 32);

}

__device__ void mont_mulModR(const uint64_t* a, const uint64_t* b, uint64_t* product)
{
  uint64_t high = 0;

  // a0 * b0
  uint128_t t = (uint128_t)a[0] * b[0] + 0;
  product[0] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a[0] * b[1] + 0 + high;
  product[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b2
  t = (uint128_t)a[0] * b[2] + 0 + high;
  product[2] = (uint64_t)t;

  // a1 * b0
  t = (uint128_t)a[1] * b[0] + product[1];
  product[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a[1] * b[1] + product[2] + high;
  product[2] = (uint64_t)t;

  // a2 * b0
  t = (uint128_t)a[2] * b[0] + product[2];
  product[2] = (uint64_t)t;


  // Mod 2^160
  product[2] &= 0xffffffff;
}

__device__ void mont_reduce(const uint64_t* t, uint64_t* y)
{
  // m1 = t_lo * k mod R
  // m2 = m1 * p / r
  // m3 = m2 + t_hi + 1

  // low 160 bits since we need to mod 2^160
  uint64_t t_lo[3];
  for(int i = 0; i < 3; i++) {
    t_lo[i] = t[i];
  }
  t_lo[2] &= 0xffffffff;

  // Remaining high bits (102 bits)
  uint64_t t_hi[3];
  t_hi[0] = (t[2] >> 32) | ((t[3] & 0xffffffff) << 32);
  t_hi[1] = (t[3] >> 32) | ((t[4] & 0xffffffff) << 32);
  t_hi[2] = 0;
  
  uint64_t m1[3];
  uint64_t m2[3];
  uint64_t m3[3];
  mont_mulModR(t_lo, _k.v, m1);

  mont_mulDivR(m1, _p.v, m2);
 
  int carry = 1;
  uint64_t sum = t_hi[0] + m2[0] + carry;
  m3[0] = sum;
  carry = sum < t_hi[0] ? 1 : 0;
  
  sum = t_hi[1] + m2[1] + carry;
  m3[1] = sum;
  carry = sum < t_hi[1] ? 1 : 0;

  sum = t_hi[2] + m2[2] + carry;
  m3[2] = sum;

  bool gte = true;

  for(int i = 2; i >= 0; i--) {
    if(m3[i] > _p.v[i]) {
      break;
    } else if(m3[i] < _p.v[i]) {
      gte = false;
      break;
    }
  }

  // Subtract P
  if(gte) {
    int borrow = 0;

    uint64_t diff = m3[0] - _p.v[0];
    borrow = diff > m3[0] ? 1 : 0;
    m3[0] = diff;
    
    diff = m3[1] - _p.v[1] - borrow;
    borrow = diff > m3[1] ? 1 : 0;
    m3[1] = diff;
    
    diff = m3[2] - _p.v[2] - borrow;
    m3[2] = diff;
  }

  for(int i = 0; i < 3; i++) {
    y[i] = m3[i];
  }
}



// Montgomery multiplication
__device__ uint131_t mul(const uint131_t x, const uint131_t y)
{
  uint64_t product[5];
  uint131_t z;
  mul160x(x.v, y.v, product);
  mont_reduce(product, z.v);

  return z;
}

__device__ uint131_t square(const uint131_t& x)
{
  uint64_t product[5];
  uint131_t z;
  square160(x.v, product);
  mont_reduce(product, z.v);

  return z;
}

__device__ uint131_t square(uint131_t x, int n)
{
  for(int i = 0; i < n; i++) {
    x = square(x);
  }

  return x;
}

__device__ bool equal(const uint131_t& x, const uint131_t& y)
{
  for(int i = 0; i < 3; i++) {
    if(x.v[i] != y.v[i]) {
      return false;
    }
  }

  return true;
}

// Modular inverse using Fermat's method
__device__ uint131_t inv_p131(uint131_t& x)
{
  uint131_t z, t0, t1, t2, t3, t4, t5;

  t0 = square(x);
  t3 = mul(x, t0);
  t4 = mul(x, t3);
  t1 = mul(x, t4);
  t2 = mul(t0, t1);
  z = mul(t0, t2);
  t4 = mul(t4, z);
  t0 = mul(t0, t4);
  t5 = mul(t3, t0);
  t5 = square(t5, 5);
  t5 = mul(t2, t5);
  t5 = square(t5, 7);
  t5 = mul(t2, t5);
  t5 = square(t5, 4);
  t5 = mul(t1, t5);
  t5 = square(t5, 8);
  t5 = mul(t0, t5);
  t5 = square(t5, 2);
  t5 = mul(t3, t5);
  t5 = square(t5, 5);
  t5 = mul(t1, t5);
  t5 = square(t5, 6);
  t5 = mul(z, t5);
  t5 = square(t5, 3);
  t5 = mul(t1, t5);
  t5 = square(t5, 7);
  t5 = mul(t4, t5);
  t5 = square(t5, 6);
  t5 = mul(t0, t5);
  t5 = square(t5, 5);
  t5 = mul(t3, t5);
  t5 = square(t5, 4);
  t5 = mul(t3, t5);
  t5 = square(t5, 5);
  t5 = mul(t3, t5);
  t5 = square(t5, 5);
  t5 = mul(t1, t5);
  t5 = square(t5, 4);
  t5 = mul(t3, t5);
  t5 = square(t5, 4);
  t5 = mul(x, t5);
  t5 = square(t5, 6);
  t5 = mul(t3, t5);
  t5 = square(t5, 5);
  t5 = mul(t3, t5);
  t5 = square(t5, 8);
  t4 = mul(t4, t5);
  t4 = square(t4, 3);
  t3 = mul(t3, t4);
  t3 = square(t3, 5);
  t2 = mul(t2, t3);
  t2 = square(t2, 4);
  t1 = mul(t1, t2);
  t1 = square(t1, 5);
  t0 = mul(t0, t1);
  t0 = square(t0, 10);
  z = mul(z, t0);

  return z;
}

__device__ uint131_t inv_p79(uint131_t& x)
{
  uint131_t z, t0, t1, t2, t3, t4, t5, t6, t7;
  
  t6 = square(x);
  t0 = mul(x, t6);
  z = square(t0);
  t3 = mul(t6, z);
  t2 = mul(t0, t3);
  t7 = mul(x, t2);
  t0 = mul(t3, t2);
  t1 = mul(t6, t0);
  t5 = mul(t6, t1);
  t4 = mul(z, t5);
  t3 = mul(t3, t4);
  t7 = mul(t7, t3);
  t6 = mul(t6, t7);
  z = mul(z, t6);
  t7 = square(t7, 7);
  t6 = mul(t6, t7);
  t6 = square(t6, 6);
  t6 = mul(t3, t6);
  t6 = square(t6, 8);
  t5 = mul(t5, t6);
  t5 = square(t5, 6);
  t4 = mul(t4, t5);
  t4 = square(t4, 11);
  t3 = mul(t3, t4);
  t3 = square(t3, 5);
  t2 = mul(t2, t3);
  t2 = square(t2, 7);
  t1 = mul(t1, t2);
  t1 = square(t1, 8);
  t0 = mul(t0, t1);
  t0 = square(t0, 8);
  t0 = mul(z, t0);
  t0 = square(t0, 6);
  z = mul(z, t0);
  z = square(z);
  z = mul(x, z);

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

__device__ bool point_exists(uint131_t& x, uint131_t& y)
{
  uint131_t y2 = square(y);
  uint131_t x3 = mul(x, square(x));
  uint131_t ax = mul(_a, x);

  uint131_t rs = add(add(x3, ax), _b);

  return equal(y2, rs);
}


#endif