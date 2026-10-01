
#ifndef _EC_MATH_CUH
#define _EC_MATH_CUH

#include "math_common.cuh"
#include "p109.cuh"
#include "p131.cuh"
#include "p79.cuh"
#include "p89.cuh"

template <int CURVE> __device__ uint131_t sub(uint131_t x, uint131_t y)
{
  return Curve<CURVE>::sub(x, y);
}

template <int CURVE> __device__ uint131_t add(const uint131_t& x, const uint131_t& y)
{
  return Curve<CURVE>::add(x, y);
}

template <int CURVE> __device__ uint131_t mul(uint131_t x, uint131_t y)
{
  return Curve<CURVE>::mul(x, y);
}

template <int CURVE> __device__ uint131_t square(uint131_t x)
{
  if constexpr(CURVE == CURVE_ID_ECP131) {
    return Curve<CURVE>::square(x);
  } else {
    return Curve<CURVE>::mul(x, x);
  }
}

template <int CURVE> __device__ uint131_t square(uint131_t x, int n)
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

  t0 = square<CURVE_ID_ECP131>(x);
  t3 = mul<CURVE_ID_ECP131>(x, t0);
  t4 = mul<CURVE_ID_ECP131>(x, t3);
  t1 = mul<CURVE_ID_ECP131>(x, t4);
  t2 = mul<CURVE_ID_ECP131>(t0, t1);
  z = mul<CURVE_ID_ECP131>(t0, t2);
  t4 = mul<CURVE_ID_ECP131>(t4, z);
  t0 = mul<CURVE_ID_ECP131>(t0, t4);
  t5 = mul<CURVE_ID_ECP131>(t3, t0);
  t5 = square<CURVE_ID_ECP131>(t5, 5);
  t5 = mul<CURVE_ID_ECP131>(t2, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 7);
  t5 = mul<CURVE_ID_ECP131>(t2, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 4);
  t5 = mul<CURVE_ID_ECP131>(t1, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 8);
  t5 = mul<CURVE_ID_ECP131>(t0, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 2);
  t5 = mul<CURVE_ID_ECP131>(t3, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 5);
  t5 = mul<CURVE_ID_ECP131>(t1, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 6);
  t5 = mul<CURVE_ID_ECP131>(z, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 3);
  t5 = mul<CURVE_ID_ECP131>(t1, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 7);
  t5 = mul<CURVE_ID_ECP131>(t4, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 6);
  t5 = mul<CURVE_ID_ECP131>(t0, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 5);
  t5 = mul<CURVE_ID_ECP131>(t3, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 4);
  t5 = mul<CURVE_ID_ECP131>(t3, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 5);
  t5 = mul<CURVE_ID_ECP131>(t3, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 5);
  t5 = mul<CURVE_ID_ECP131>(t1, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 4);
  t5 = mul<CURVE_ID_ECP131>(t3, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 4);
  t5 = mul<CURVE_ID_ECP131>(x, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 6);
  t5 = mul<CURVE_ID_ECP131>(t3, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 5);
  t5 = mul<CURVE_ID_ECP131>(t3, t5);
  t5 = square<CURVE_ID_ECP131>(t5, 8);
  t4 = mul<CURVE_ID_ECP131>(t4, t5);
  t4 = square<CURVE_ID_ECP131>(t4, 3);
  t3 = mul<CURVE_ID_ECP131>(t3, t4);
  t3 = square<CURVE_ID_ECP131>(t3, 5);
  t2 = mul<CURVE_ID_ECP131>(t2, t3);
  t2 = square<CURVE_ID_ECP131>(t2, 4);
  t1 = mul<CURVE_ID_ECP131>(t1, t2);
  t1 = square<CURVE_ID_ECP131>(t1, 5);
  t0 = mul<CURVE_ID_ECP131>(t0, t1);
  t0 = square<CURVE_ID_ECP131>(t0, 10);
  z = mul<CURVE_ID_ECP131>(z, t0);

  return z;
}

__device__ uint131_t inv_p79(uint131_t& x)
{
  uint131_t z, t0, t1, t2, t3, t4, t5, t6, t7;

  t6 = square<CURVE_ID_ECP79>(x);
  t0 = mul<CURVE_ID_ECP79>(x, t6);
  z = square<CURVE_ID_ECP79>(t0);
  t3 = mul<CURVE_ID_ECP79>(t6, z);
  t2 = mul<CURVE_ID_ECP79>(t0, t3);
  t7 = mul<CURVE_ID_ECP79>(x, t2);
  t0 = mul<CURVE_ID_ECP79>(t3, t2);
  t1 = mul<CURVE_ID_ECP79>(t6, t0);
  t5 = mul<CURVE_ID_ECP79>(t6, t1);
  t4 = mul<CURVE_ID_ECP79>(z, t5);
  t3 = mul<CURVE_ID_ECP79>(t3, t4);
  t7 = mul<CURVE_ID_ECP79>(t7, t3);
  t6 = mul<CURVE_ID_ECP79>(t6, t7);
  z = mul<CURVE_ID_ECP79>(z, t6);
  t7 = square<CURVE_ID_ECP79>(t7, 7);
  t6 = mul<CURVE_ID_ECP79>(t6, t7);
  t6 = square<CURVE_ID_ECP79>(t6, 6);
  t6 = mul<CURVE_ID_ECP79>(t3, t6);
  t6 = square<CURVE_ID_ECP79>(t6, 8);
  t5 = mul<CURVE_ID_ECP79>(t5, t6);
  t5 = square<CURVE_ID_ECP79>(t5, 6);
  t4 = mul<CURVE_ID_ECP79>(t4, t5);
  t4 = square<CURVE_ID_ECP79>(t4, 11);
  t3 = mul<CURVE_ID_ECP79>(t3, t4);
  t3 = square<CURVE_ID_ECP79>(t3, 5);
  t2 = mul<CURVE_ID_ECP79>(t2, t3);
  t2 = square<CURVE_ID_ECP79>(t2, 7);
  t1 = mul<CURVE_ID_ECP79>(t1, t2);
  t1 = square<CURVE_ID_ECP79>(t1, 8);
  t0 = mul<CURVE_ID_ECP79>(t0, t1);
  t0 = square<CURVE_ID_ECP79>(t0, 8);
  t0 = mul<CURVE_ID_ECP79>(z, t0);
  t0 = square<CURVE_ID_ECP79>(t0, 6);
  z = mul<CURVE_ID_ECP79>(z, t0);
  z = square<CURVE_ID_ECP79>(z);
  z = mul<CURVE_ID_ECP79>(x, z);

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
      prod = mul<CURVE_ID_ECP89>(prod, y);
    }
    y = square<CURVE_ID_ECP89>(y);

    bits >>= 1;
  }

  bits = _p89_p.w.v1;
  for(int i = 0; i < 25; i++) {
    if(bits & 1) {
      prod = mul<CURVE_ID_ECP89>(prod, y);
    }
    y = square<CURVE_ID_ECP89>(y);

    bits >>= 1;
  }

  return prod;
}

__device__ uint131_t inv_p109(uint131_t& x)
{
  uint131_t t10 = square<CURVE_ID_ECP109>(x);
  uint131_t t11 = mul<CURVE_ID_ECP109>(x, t10);
  uint131_t t101 = mul<CURVE_ID_ECP109>(t10, t11);
  uint131_t t111 = mul<CURVE_ID_ECP109>(t10, t101);
  uint131_t t1001 = mul<CURVE_ID_ECP109>(t10, t111);
  uint131_t t1011 = mul<CURVE_ID_ECP109>(t10, t1001);
  uint131_t t1101 = mul<CURVE_ID_ECP109>(t10, t1011);
  uint131_t t1111 = mul<CURVE_ID_ECP109>(t10, t1101);
  uint131_t t11010 = mul<CURVE_ID_ECP109>(t1011, t1111);
  uint131_t t1101000 = square<CURVE_ID_ECP109>(t11010, 2);
  uint131_t t1101111 = mul<CURVE_ID_ECP109>(t111, t1101000);

  uint131_t i25 = square<CURVE_ID_ECP109>(t1101111, 4);
  i25 = mul<CURVE_ID_ECP109>(i25, t101);
  i25 = square<CURVE_ID_ECP109>(i25, 5);
  i25 = mul<CURVE_ID_ECP109>(i25, t1011);
  i25 = square<CURVE_ID_ECP109>(i25, 2);

  uint131_t i36 = mul<CURVE_ID_ECP109>(t11, i25);
  i36 = square<CURVE_ID_ECP109>(i36, 6);
  i36 = mul<CURVE_ID_ECP109>(i36, t1011);
  i36 = square<CURVE_ID_ECP109>(i36, 2);
  i36 = mul<CURVE_ID_ECP109>(i36, t11);

  uint131_t i54 = square<CURVE_ID_ECP109>(i36, 6);
  i54 = mul<CURVE_ID_ECP109>(i54, t1001);
  i54 = square<CURVE_ID_ECP109>(i54, 5);
  i54 = mul<CURVE_ID_ECP109>(i54, t1011);
  i54 = square<CURVE_ID_ECP109>(i54, 5);

  uint131_t i73 = mul<CURVE_ID_ECP109>(t111, i54);
  i73 = square<CURVE_ID_ECP109>(i73, 11);
  i73 = mul<CURVE_ID_ECP109>(i73, t1011);
  i73 = square<CURVE_ID_ECP109>(i73, 5);
  i73 = mul<CURVE_ID_ECP109>(i73, t1011);

  uint131_t i93 = square<CURVE_ID_ECP109>(i73, 5);
  i93 = mul<CURVE_ID_ECP109>(i93, t1101);
  i93 = square<CURVE_ID_ECP109>(i93, 5);
  i93 = mul<CURVE_ID_ECP109>(i93, t1001);
  i93 = square<CURVE_ID_ECP109>(i93, 8);

  uint131_t i106 = mul<CURVE_ID_ECP109>(t1111, i93);
  i106 = square<CURVE_ID_ECP109>(i106, 6);
  i106 = mul<CURVE_ID_ECP109>(i106, t1101);
  i106 = square<CURVE_ID_ECP109>(i106, 4);
  i106 = mul<CURVE_ID_ECP109>(i106, t1011);

  uint131_t i121 = square<CURVE_ID_ECP109>(i106, 6);
  i121 = mul<CURVE_ID_ECP109>(i121, t1111);
  i121 = square<CURVE_ID_ECP109>(i121, 4);
  i121 = mul<CURVE_ID_ECP109>(i121, t1101);
  i121 = square<CURVE_ID_ECP109>(i121, 3);

  uint131_t i133 = mul<CURVE_ID_ECP109>(t101, i121);
  i133 = square<CURVE_ID_ECP109>(i133, 3);
  i133 = mul<CURVE_ID_ECP109>(i133, t11);
  i133 = square<CURVE_ID_ECP109>(i133, 6);
  i133 = mul<CURVE_ID_ECP109>(i133, t1011);

  return mul<CURVE_ID_ECP109>(square<CURVE_ID_ECP109>(i133), x);
}

template <int CURVE> __device__ uint131_t inv(uint131_t x)
{
  uint131_t r;
  if constexpr(CURVE == CURVE_ID_ECP131) {
    r = inv_p131(x);
  } else if constexpr(CURVE == CURVE_ID_ECP79) {
    r = inv_p79(x);
  } else if constexpr(CURVE == CURVE_ID_ECP89) {
    r = inv_p89(x);
  } else if constexpr(CURVE == CURVE_ID_ECP109) {
    r = inv_p109(x);
  }

  return r;
}

// ECC functons

// Check for point-at-infinity using the x coordinate
__device__ bool is_infinity(uint131_t x) { return (x.w.v2 & 0xff) == 0xff; }

__device__ void set_point_at_infinity(uint131_t& x) { x.w.v2 = (uint32_t)-1; }

template <int CURVE> __device__ bool point_exists_prime(uint131_t& x, uint131_t& y)
{
  uint131_t a = Curve<CURVE>::a();
  uint131_t b = Curve<CURVE>::b();

  uint131_t y2 = square<CURVE>(y);
  uint131_t x3 = mul<CURVE>(x, square<CURVE>(x));
  uint131_t ax = mul<CURVE>(a, x);

  uint131_t rs = add<CURVE>(add<CURVE>(x3, ax), b);

  return equal(y2, rs);
}

template <int CURVE> __device__ bool point_exists(uint131_t& x, uint131_t& y)
{
  return point_exists_prime<CURVE>(x, y);
}

#endif
