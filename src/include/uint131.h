#ifndef _UINT131_H
#define _UINT131_H

#include <stdint.h>
#include <string.h>
#include <string>

#include "shared_types.h"

inline uint131_t sub_raw(const uint131_t& x, const uint131_t& y)
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

inline uint131_t add_raw(const uint131_t& x, const uint131_t& y, int carry_in = 0)
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

inline uint131_t rshift(uint131_t x, int n)
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
      y.v[i] |= (x.v[i + off + 1] << (left_shift));
    }
  }

  return y;
}

inline uint131_t lshift(uint131_t x, int n)
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

// if x is less than y
inline int is_less_than(uint131_t& x, uint131_t& y)
{
  uint131_t diff = sub_raw(x, y);

  return (diff.w.v2 >> 31) & 1;
}

uint131_t make_uint131(uint32_t x);
uint131_t make_uint131(const std::string& hex);

// bool operator==(const uint131_t& a, const uint131_t& b);
// bool operator!=(const uint131_t& a, const uint131_t& b);
// bool is_odd(const uint131_t& x);

inline bool operator==(const uint131_t& a, const uint131_t& b)
{
  return a.w.v0 == b.w.v0 && a.w.v1 == b.w.v1 && a.w.v2 == b.w.v2;
}

inline bool operator!=(const uint131_t& a, const uint131_t& b) { return !(a == b); }

inline bool is_odd(const uint131_t& x) { return x.w.v0 & 0x01; }

inline uint131_t add_mod_n(uint131_t x, uint131_t y, uint131_t n)
{
  uint131_t z = add_raw(x, y);

  // Reduce mod P
  if(is_less_than(n, z)) {
    z = sub_raw(z, n);
  }
  return z;
}

inline uint131_t sub_mod_n(uint131_t x, uint131_t y, uint131_t n)
{
  uint131_t z = sub_raw(x, y);
  int borrow = z.w.v2 >> 31;

  // Went below zero. Need to add P.
  if(borrow) {
    z = add_raw(z, n);
  }

  return z;
}

inline uint131_t mul_mod_n(uint131_t x, uint131_t y, uint131_t n)
{
  uint131_t product = make_uint131(0);

  for(int i = 0; i < 131; i++) {
    int word = i / 32;
    int bit = i % 32;

    if(x.v[word] & (1 << bit)) {
      product = add_raw(product, y);
    }
    if(is_less_than(n, product)) {
      product = sub_raw(product, n);
    }

    y = lshift(y, 1);
    if(is_less_than(n, y)) {
      y = sub_raw(y, n);
    }
  }

  return product;
}

inline uint131_t inv_mod_n(uint131_t x, uint131_t n)
{
  uint131_t m = sub_raw(n, make_uint131(2));
  uint131_t product = make_uint131(1);

  for(int i = 0; i < 131; i++) {
    if(m.v[0] & 1) {
      product = mul_mod_n(product, x, n);
    }
    m = rshift(m, 1);
    x = mul_mod_n(x, x, n);
  }

  return product;
}

std::string to_str(const uint131_t& x);

#endif