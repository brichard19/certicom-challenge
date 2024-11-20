#include <stdexcept>

#include "eccp131.h"
#include "mont.h"

typedef unsigned __int128 uint128_t;

uint131_t _p;
uint131_t _a;
uint131_t _b;
uint131_t _n;
uint131_t _gx;
uint131_t _gy;
uint131_t _qx;
uint131_t _qy;
uint131_t _k;
uint131_t _one;
uint131_t _p_minus_2;
uint131_t _r;
uint131_t _r2;
uint131_t _modulo_sqrt;

int _bits;
int _words;
std::string _curve_name;

namespace ecc {
void set_curve(const std::string& curve_name)
{
  if(curve_name == "ecp131") {
    _p = _eccp131_p;
    _a = _eccp131_a;
    _b = _eccp131_b;
    _n = _eccp131_n;
    _gx = _eccp131_gx;
    _gy = _eccp131_gy;
    _qx = _eccp131_qx;
    _qy = _eccp131_qy;
    _k = _eccp131_k;
    _modulo_sqrt = _eccp131_sqrt;
    _one = _eccp131_one;
    _p_minus_2 = _eccp131_p_minus_2;
    _r = _eccp131_r;
    _r2 = _eccp131_r2;
    _bits = _eccp131_bits;
    _words = (_bits + 63) / 64;
    _curve_name = _eccp131_curve_name;
  } else if(curve_name == "ecp79") {
    _p = _eccp79_p;
    _a = _eccp79_a;
    _b = _eccp79_b;
    _n = _eccp79_n;
    _gx = _eccp79_gx;
    _gy = _eccp79_gy;
    _qx = _eccp79_qx;
    _qy = _eccp79_qy;
    _k = _eccp79_k;
    _modulo_sqrt = _eccp79_sqrt;
    _one = _eccp79_one;
    _p_minus_2 = _eccp79_p_minus_2;
    _r = _eccp79_r;
    _r2 = _eccp79_r2;
    _bits = _eccp79_bits;
    _words = (_bits + 63) / 64;
    _curve_name = _eccp79_curve_name;
  } else {
    throw std::runtime_error("Invalid curve name");
  }
}

}

namespace {

// 131 x 131 -> 262 multiplication
void mul_131(const uint64_t* a, const uint64_t* b, uint64_t* product)
{
  uint64_t tmp[5] = {0};
  uint64_t high = 0;

  // a0 * b0
  uint128_t t = (uint128_t)a[0] * b[0] + tmp[0];
  tmp[0] = (uint64_t)t;
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

  product[0] = tmp[0];
  product[1] = tmp[1];
  product[2] = tmp[2];
  product[3] = tmp[3];
  product[4] = tmp[4];
}

// 131 x 131 -> 262 multiplication
void square_131(const uint64_t* a, uint64_t* product)
{
  uint64_t tmp[5] = {0};
  uint64_t high = 0;

  // a0 * a0
  uint128_t t = (uint128_t)a[0] * a[0] + tmp[0];
  tmp[0] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * a1 
  t = (uint128_t)a[0] * a[1] + tmp[1] + high;
  tmp[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * a2
  t = (uint128_t)a[0] * a[2] + tmp[2] + high;
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  tmp[3] = high;

  // a1 * a0
  t = (uint128_t)a[1] * a[0] + tmp[1];
  tmp[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * a1 
  t = (uint128_t)a[1] * a[1] + tmp[2] + high;
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a1 * a2
  t = (uint128_t)a[1] * a[2] + tmp[3] + high;
  tmp[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  tmp[4] = high;

  // a2 * a0
  t = (uint128_t)a[2] * a[0] + tmp[2];
  tmp[2] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * a1 
  t = (uint128_t)a[2] * a[1] + tmp[3] + high;
  tmp[3] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a2 * a2
  // The final word is only at most 6 bits, so no 128-bit mul needed
  uint64_t t64 = a[2] * a[2] + tmp[4] + high;
  tmp[4] = t64;

  product[0] = tmp[0];
  product[1] = tmp[1];
  product[2] = tmp[2];
  product[3] = tmp[3];
  product[4] = tmp[4];
}

// One in put is 131 bits, the other is 160 bits
// Output is 131 bits
void mul_shift_160(const uint64_t* a, const uint64_t* b, uint64_t* product)
{
  uint64_t tmp[5] = {0};
  uint64_t high = 0;

  // a0 * b0
  uint128_t t = (uint128_t)a[0] * b[0] + tmp[0];
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
  // The final word is only at most 35 bits, so no 128-bit mul needed
  uint64_t t64 = a[2] * b[2] + tmp[4] + high;
  tmp[4] = t64;

  // Divide by 2^160
  product[0] = (tmp[2] >> 32) | ((tmp[3] & 0xffffffff) << 32);
  product[1] = (tmp[3] >> 32) | ((tmp[4] & 0xffffffff) << 32);
  product[2] = (tmp[4] >> 32);
}

void mul_mod_160(const uint64_t* a, const uint64_t* b, uint64_t* product)
{
  uint64_t tmp[3] = {0};
  uint64_t high = 0;

  // a0 * b0
  uint128_t t = (uint128_t)a[0] * b[0] + tmp[0];
  tmp[0] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a[0] * b[1] + tmp[1] + high;
  tmp[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b2
  t = (uint128_t)a[0] * b[2] + tmp[2] + high;
  tmp[2] = (uint64_t)t;

  // a1 * b0
  t = (uint128_t)a[1] * b[0] + tmp[1];
  tmp[1] = (uint64_t)t;
  high = (uint64_t)(t >> 64);

  // a0 * b1 
  t = (uint128_t)a[1] * b[1] + tmp[2] + high;
  tmp[2] = (uint64_t)t;

  // a2 * b0
  t = (uint128_t)a[2] * b[0] + tmp[2];
  tmp[2] = (uint64_t)t;

  // Mod 2^160
  product[0] = tmp[0];
  product[1] = tmp[1];
  product[2] = tmp[2] & 0xffffffff;
}

void mont_reduce(const uint64_t* t, uint64_t* y)
{
  // m1 = t_lo * k mod R
  // m2 = m1 * p / r
  // m3 = m2 + t_hi + 1

  // low 160 bits since we need to do mod 2^160
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
  mul_mod_160(t_lo, _k.v, m1);

  mul_shift_160(m1, _p.v, m2);
 
  // Not sure if this is correct, but carry always seems to be 1.
  uint64_t carry = 1;

  for(int i = 0; i < 3; i++) {
    uint128_t sum = (uint128_t)t_hi[i] + m2[i] + carry;

    m3[i] = (uint64_t)sum;
    carry = (uint64_t)(sum >> 64);
  }

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
    uint64_t borrow = 0;
    for(int i = 0; i < 3; i++) {
      uint128_t diff = (uint128_t)m3[i] - _p.v[i] - borrow;
      borrow = (uint64_t)(diff >> 64) & 1;
      m3[i] = (uint64_t)diff;
    }
  }

  for(int i = 0; i < 3; i++) {
    y[i] = m3[i];
  }
}

}


namespace mont {


uint131_t to(uint131_t x)
{
  return mul(x, _r2);
}

uint131_t from(uint131_t x)
{
  uint131_t one = {{1, 0, 0}};
  return mul(x, one);
}



uint131_t add(uint131_t x, uint131_t y)
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

uint131_t neg(uint131_t x)
{
  // To get -x mod P, subtract x from P

  uint131_t z = {{0}};

  uint64_t borrow = 0;

  for(int i = 0; i < 3; i++) {
    uint128_t diff = (uint128_t)x.v[i] - _p.v[i] - borrow;

    z.v[i] = (uint64_t)diff;
    borrow = (uint64_t)(diff >> 64) & 1;
  }

  return z;
}

uint131_t sub(uint131_t x, uint131_t y)
{
  uint131_t z = {{0}};

  uint64_t borrow = 0;

  for(int i = 0; i < 3; i++) {
    uint128_t diff = (uint128_t)x.v[i] - y.v[i] - borrow;

    z.v[i] = (uint64_t)diff;
    borrow = (uint64_t)(diff >> 64) & 1;
  }

  // Went below zero. Need to add P.
  if(borrow) {
    uint64_t carry = 0;
    for(int i = 0; i < 3; i++) {
      uint128_t sum = (uint128_t)z.v[i] + _p.v[i] + carry;

      z.v[i] = (uint64_t)sum;
      carry = (uint64_t)(sum >> 64);
    }
  }

  return z;
}

uint131_t mul(uint131_t x, uint131_t y)
{
  uint64_t product[5] = {0};
  uint131_t z;

  mul_131(x.v, y.v, product);
  mont_reduce(product, z.v);

  return z;
}

uint131_t square(uint131_t x)
{
  uint64_t product[5] = {0};
  uint131_t z;

  square_131(x.v, product);
  mont_reduce(product, z.v);

  return z;
}

uint131_t pow(uint131_t x, uint131_t exponent)
{
    // Initialize to 1 in montgomery form
    uint131_t product = _one;
    uint131_t q = x;

    for(int w = 0; w < 3; w++) {
      uint64_t mask = 0x00000001;

      while(mask) {
        if(exponent.v[w] & mask) {
          product = mul(product, q);
        }
        q = square(q);

        mask <<= 1;
      }
    }

    return product;
}

// Modular inverse using Fermat's method
uint131_t inv(uint131_t x)
{
  // Initialize to 1 in montgomery form
  uint131_t product = _one;
  uint131_t q = x;

  for(int w = 0; w < 3; w++) {
    uint64_t mask = 0x00000001;

    while(mask) {
      if(_p_minus_2.v[w] & mask) {
        product = mul(product, q);
      }
      q = square(q);

      mask <<= 1;
    }
  }

  return product;
}

uint131_t sqrt(uint131_t x)
{
  if(_curve_name == "ecp131") {

    // sqrt(x) = +- x^(p+1)/4 for p congruent to 3 mod 4
    return pow(x, _modulo_sqrt);
  } else if(_curve_name == "ecp79") {

    // 2 in montgomery form 
    uint131_t two = {{0xa88f54e07ed578be, 0x26b0, 0x00}};
    
    uint131_t v = pow(mul(two, x), _modulo_sqrt);

    uint131_t i = mul(mul(square(v), x), two);

    uint131_t xv = mul(x, v);

    return mul(xv, sub(i, _one));
  } else {
    throw std::runtime_error("Invalid curve name");
  }
}

}