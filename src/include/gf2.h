#ifndef GF2_H
#define GF2_H

#include "uint131.h"

namespace gf2 {

// Addition and subtraction are both coefficient-wise XOR in GF(2).
inline uint131_t add(const uint131_t& a, const uint131_t& b)
{
  uint131_t result = {};
  for(int i = 0; i < 5; i++) {
    result.v[i] = a.v[i] ^ b.v[i];
  }
  return result;
}

namespace detail {

inline bool bit(const uint131_t& value, int index)
{
  return (value.v[index / 32] & (uint32_t(1) << (index % 32))) != 0;
}

// Returns the degree of a non-zero polynomial, or -1 for zero.
inline int degree(const uint131_t& polynomial)
{
  for(int i = 159; i >= 0; i--) {
    if(bit(polynomial, i)) {
      return i;
    }
  }
  return -1;
}

} // namespace detail

// Multiplies two elements modulo a polynomial of degree at most 131.
// The modulus includes its leading x^n term, e.g. x^3 + x + 1 is 0b1011.
// A zero or constant modulus is invalid and returns zero.
inline uint131_t mul(uint131_t a, uint131_t b, const uint131_t& modulus)
{
  const int n = detail::degree(modulus);
  if(n < 1 || n > 131) {
    return {};
  }

  uint131_t product = {};

  // Reduce the operands first so callers may pass any polynomials that fit
  // uint131_t, rather than only canonical field elements.
  for(int i = 159; i >= n; i--) {
    if(detail::bit(a, i)) {
      a = add(a, lshift(modulus, i - n));
    }
    if(detail::bit(b, i)) {
      b = add(b, lshift(modulus, i - n));
    }
  }

  for(int i = 0; i < n; i++) {
    if(detail::bit(b, i)) {
      product = add(product, a);
    }

    const bool reduce = detail::bit(a, n - 1);
    a = lshift(a, 1);
    if(reduce) {
      a = add(a, modulus);
    }
  }

  return product;
}

// Returns the multiplicative inverse of a modulo an irreducible polynomial.
// Zero has no inverse and returns zero.
inline uint131_t inv(const uint131_t& a, const uint131_t& modulus)
{
  const int n = detail::degree(modulus);
  if(n < 1 || n > 131 || detail::degree(a) < 0) {
    return {};
  }

  // In GF(2^n), a^(2^n - 1) = 1 for non-zero a, so the inverse is
  // a^(2^n - 2). The exponent has n - 1 one bits followed by a zero bit.
  uint131_t result = make_uint131(1);
  for(int i = n - 1; i >= 0; i--) {
    result = mul(result, result, modulus);
    if(i != 0) {
      result = mul(result, a, modulus);
    }
  }

  return result;
}

} // namespace gf2

#endif
