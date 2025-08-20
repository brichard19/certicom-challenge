#include <stdexcept>

#include "eccp131.h"
#include "montgomery.h"

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

uint131_t sub_raw(const uint131_t& x, const uint131_t& y)
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

  return z;
}

// if x is less than y
int is_less_than(uint131_t& x, uint131_t& y)
{
  uint131_t diff = sub_raw(x, y);

  return (diff.w.v2 >> 31) & 1;
}

uint131_t add_raw(const uint131_t& x, const uint131_t& y, int carry_in = 0)
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

  return z;
}

// 131 x 131 -> 262 multiplication
uint262_t mul_131(const uint131_t& a, const uint131_t& b)
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
uint262_t square_131(const uint131_t& a)
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
uint131_t mul_shift_160(const uint160_t& a, const uint131_t& b)
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

//void mul_mod_160(const uint64_t* a, const uint64_t* b, uint64_t* product)
uint160_t mul_mod_160(const uint160_t& a, const uint131_t& b)
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

//void mont_reduce(const uint64_t* t, uint64_t* y)
uint131_t mont_reduce(const uint262_t& t)
{
  // m1 = t_lo * k mod R
  // m2 = m1 * p / r
  // m3 = m2 + t_hi + 1

  uint160_t t_lo; 
  t_lo.w.v0 = t.v[0];
  t_lo.w.v1 = t.v[1];
  t_lo.w.v2 = (uint32_t)t.v[2];

  // Remaining high bits (102 bits)
  uint131_t t_hi;

  t_hi.w.v0 = (t.v[2] >> 32) | ((t.v[3] & 0xffffffff) << 32);
  t_hi.w.v1 = (t.v[3] >> 32) | ((t.v[4] & 0xffffffff) << 32);
  t_hi.w.v2 = 0;

  uint160_t m1;
  uint131_t m2;
  uint131_t m3;

  m1 = mul_mod_160(t_lo, _k);

  m2 = mul_shift_160(m1, _p);
 
  // Not sure if this is correct, but carry always seems to be 1.
  m3 = add_raw(t_hi, m2, 1);

  if(is_less_than(_p, m3)) {
    m3 = sub_raw(m3, _p);
  }
  return m3;

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
  uint131_t z = add_raw(x, y);

  // Reduce mod P
  if(is_less_than(_p, z)) {
    z = sub_raw(z, _p);
  }
  return z;
}

uint131_t neg(uint131_t x)
{
  // To get -x mod P, subtract x from P
  return sub_raw(_p, x);
}



uint131_t sub(uint131_t x, uint131_t y)
{
  uint131_t z = sub_raw(x, y);
  int borrow = z.w.v2 >> 31;

  // Went below zero. Need to add P.
  if(borrow) {
    z = add_raw(z, _p);
  }

  return z;
}

uint131_t rshift(uint131_t x, int n)
{
  if(n == 0) {
    return x;
  } 

  uint131_t y = {0};
  int off = n / 32;
  int right_shift = n % 32;
  int left_shift = right_shift == 0 ? 0 : 32 - right_shift;

  for(int i = 0; i + off < 5; i++) {
    y.v[i] = (x.v[i + off] >> right_shift);
    if(right_shift != 0 && i + off + 1 < 5) {
      y.v[i] |=(x.v[i + off + 1] << (left_shift));
    }
  }

  return y;
}

uint131_t lshift(uint131_t x, int n)
{
  if(n == 0) {
    return x;
  }

  uint131_t y = {0};

  int off = n / 32;
  int left_shift = n % 32;
  int right_shift = left_shift == 0 ? 0 : 32 - left_shift;

  for(int i = 4; i - off >= 0; i--) {
    y.v[i] = (x.v[i - off] << left_shift);
    if(left_shift != 0 && i - off - 1 >= 0) {
      y.v[i] |= (x.v[i - off - 1] >> (right_shift));
    }
  }

  return y;
}


uint131_t mul(uint131_t x, uint131_t y)
{
  uint262_t product = mul_131(x, y);
  return mont_reduce(product);
}

uint131_t square(uint131_t x)
{
  uint262_t product = square_131(x);
  return mont_reduce(product);
}

uint131_t pow(uint131_t x, uint131_t exponent)
{
    // Initialize to 1 in montgomery form
    uint131_t product = _one;
    uint131_t q = x;

    for(int i = 0; i < 64; i++) {
      if(exponent.w.v0 & 1) {
        product = mul(product, q);
      }
      q = square(q);

      exponent.w.v0 >>= 1;
    }
    
    for(int i = 0; i < 64; i++) {
      if(exponent.w.v1 & 1) {
        product = mul(product, q);
      }
      q = square(q);

      exponent.w.v1 >>= 1;
    }

    while(exponent.w.v2) {
      if(exponent.w.v2 & 1) {
        product = mul(product, q);
      }
      q = square(q);
      exponent.w.v2 >>= 1;
    }

    return product;
}

// Modular inverse using Fermat's method
uint131_t inv(uint131_t x)
{
  return pow(x, _p_minus_2);
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