
#ifndef _EC_MATH_CUH
#define _EC_MATH_CUH

#include "math_common.cuh"
#include "p131.cuh"
#include "p79.cuh"
#include "p89.cuh"


template<int CURVE> __device__ uint131_t sub(uint131_t x, uint131_t y)
{
  return Curve<CURVE>::sub(x, y);
}


template<int CURVE> __device__ uint131_t add(const uint131_t& x, const uint131_t& y)
{
  return Curve<CURVE>::add(x, y);
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

  // Remaining high bits (102 bits)
  uint131_t t_hi = Curve<CURVE>::high_bits(t);

  uint160_t m1;
  uint131_t m2;
  uint131_t m3;

  uint131_t p = Curve<CURVE>::p();
  uint131_t k = Curve<CURVE>::k();

  m1 = mul_mod_160(t_lo, k);

  m2 = Curve<CURVE>::mul_shift_160(m1, p);

  // Not sure if this is correct, but carry always seems to be 1.
  m3 = add_raw(t_hi, m2, 1);

  if(is_less_than(p, m3)) {
    m3 = sub_raw(m3, p);
  }

  return m3;
}


template<int CURVE> __device__ uint131_t mul(uint131_t x, uint131_t y)
{
  uint262_t product = Curve<CURVE>::mul(x, y);
  
  return mont_reduce<CURVE>(product);
}


template<int CURVE> __device__ uint131_t square(uint131_t x)
{
  uint262_t product = Curve<CURVE>::square(x);

  return mont_reduce<CURVE>(product);
}


template<int CURVE> __device__ uint131_t square(uint131_t x, int n)
{
  for(int i = 0; i < n; i++) {
    x = square<CURVE>(x);
  }

  return x;
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

// TODO: Optimize
__device__ uint131_t inv_p89(uint131_t& x)
{
  uint131_t prod = _p89_one;
  uint131_t y = x;

  uint64_t bits = _p89_p.w.v0 - 2;

  for(int i = 0; i < 64; i++) {
    if(bits & 1) {
      prod = mul<89>(prod, y);
    }
    y = square<89>(y);
    
    bits >>= 1;
  }

  bits = _p89_p.w.v1;
  for(int i = 0; i < 25; i++) {
    if(bits & 1) {
      prod = mul<89>(prod, y);
    }
    y = square<89>(y);
    
    bits >>= 1;
  }

  return prod;
}

template<int CURVE> __device__ uint131_t inv(uint131_t x)
{
  uint131_t r;
  if(CURVE == 131) {
    r = inv_p131(x);
  } else if(CURVE == 79) {
    r = inv_p79(x);
  } else if(CURVE == 89) {
    r = inv_p89(x);
  }
  
  return r;
}

// ECC functons

// Check for point-at-infinity using the x coordinate
__device__ bool is_infinity(uint131_t x)
{
  return (x.w.v2 & 0xff) == 0xff;
}

__device__ void set_point_at_infinity(uint131_t& x)
{
  x.w.v2 = (uint32_t)-1;
}

template<int CURVE> __device__ bool point_exists(uint131_t& x, uint131_t& y)
{
  uint131_t a = Curve<CURVE>::a();
  uint131_t b = Curve<CURVE>::b();

  uint131_t y2 = square<CURVE>(y);
  uint131_t x3 = mul<CURVE>(x, square<CURVE>(x));
  uint131_t ax = mul<CURVE>(a, x);

  uint131_t rs = add<CURVE>(add<CURVE>(x3, ax), b);

  return equal(y2, rs);
}


#endif