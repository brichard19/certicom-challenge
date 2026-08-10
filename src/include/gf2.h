#ifndef GF2_H
#define GF2_H

#include "uint131.h"

#include <immintrin.h>
#include <stdexcept>

namespace gf2 {

enum class Field {
  NONE,
  GF2_79,
  GF2_89,
  GF2_131,
};

// Addition and subtraction are both coefficient-wise XOR in GF(2).
inline uint131_t add(const uint131_t& a, const uint131_t& b)
{
  uint131_t result = {};
  for(int i = 0; i < 5; i++)
    result.v[i] = a.v[i] ^ b.v[i];
  return result;
}

namespace detail {

#if defined(__GNUC__) || defined(__clang__)
__attribute__((target("pclmul,sse2")))
#endif
inline void
clmul64(uint64_t a, uint64_t b, uint64_t& low, uint64_t& high)
{
  __m128i product = _mm_clmulepi64_si128(_mm_set_epi64x(0, a), _mm_set_epi64x(0, b), 0x00);
  low = uint64_t(_mm_cvtsi128_si64(product));
  high = uint64_t(_mm_cvtsi128_si64(_mm_srli_si128(product, 8)));
}

inline void multiply(const uint131_t& a, const uint131_t& b, uint64_t product[6])
{
  const uint64_t aw[3] = {a.w.v0, a.w.v1, a.w.v2};
  const uint64_t bw[3] = {b.w.v0, b.w.v1, b.w.v2};

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      uint64_t low;
      uint64_t high;
      clmul64(aw[i], bw[j], low, high);
      product[i + j] ^= low;
      product[i + j + 1] ^= high;
    }
  }
}

inline uint131_t reduce_131(const uint64_t product[6])
{
  uint64_t h0 = (product[2] >> 3) | (product[3] << 61);
  uint64_t h1 = (product[3] >> 3) | (product[4] << 61);
  uint64_t h2 = product[4] >> 3;

  uint64_t r0 = product[0] ^ h0 ^ (h0 << 1) ^ (h0 << 2) ^ (h0 << 13);
  uint64_t r1 =
      product[1] ^ h1 ^ (h1 << 1) ^ (h0 >> 63) ^ (h1 << 2) ^ (h0 >> 62) ^ (h1 << 13) ^ (h0 >> 51);
  uint64_t r2 = (product[2] & 0x7) ^ h2 ^ (h2 << 1) ^ (h1 >> 63) ^ (h2 << 2) ^ (h1 >> 62) ^
                (h2 << 13) ^ (h1 >> 51);

  uint64_t high = r2 >> 3;
  r2 &= 0x7;
  r0 ^= high ^ (high << 1) ^ (high << 2) ^ (high << 13);

  uint131_t result = {};
  result.w.v0 = r0;
  result.w.v1 = r1;
  result.w.v2 = uint32_t(r2);
  return result;
}

inline uint131_t reduce_89(const uint64_t product[6])
{
  uint64_t h0 = (product[1] >> 25) | (product[2] << 39);
  uint64_t h1 = product[2] >> 25;

  uint64_t r0 = product[0] ^ h0 ^ (h0 << 38);
  uint64_t r1 = (product[1] & 0x1ffffff) ^ h1 ^ (h0 >> 26) ^ (h1 << 38);

  uint64_t high = r1 >> 25;
  r1 &= 0x1ffffff;
  r0 ^= high ^ (high << 38);
  r1 ^= high >> 26;

  uint131_t result = {};
  result.w.v0 = r0;
  result.w.v1 = r1;
  return result;
}

inline uint131_t reduce_79(const uint64_t product[6])
{
  uint64_t h0 = (product[1] >> 15) | (product[2] << 49);
  uint64_t h1 = product[2] >> 15;

  uint64_t r0 = product[0] ^ h0 ^ (h0 << 9);
  uint64_t r1 = (product[1] & 0x7fff) ^ h1 ^ (h0 >> 55) ^ (h1 << 9);

  uint64_t high = r1 >> 15;
  r1 &= 0x7fff;
  r0 ^= high ^ (high << 9);

  uint131_t result = {};
  result.w.v0 = r0;
  result.w.v1 = r1;
  return result;
}

} // namespace detail

inline uint131_t mul(const uint131_t& a, const uint131_t& b, Field field)
{
  uint64_t product[6] = {};
  detail::multiply(a, b, product);

  switch(field) {
  case Field::GF2_131:
    return detail::reduce_131(product);
  case Field::GF2_89:
    return detail::reduce_89(product);
  case Field::GF2_79:
    return detail::reduce_79(product);
  default:
    throw std::invalid_argument("Unsupported binary field");
  }
}

inline uint131_t inv(const uint131_t& a, Field field)
{
  if(a == make_uint131(0))
    return {};

  int degree;
  if(field == Field::GF2_131) {
    degree = 131;
  } else if(field == Field::GF2_89) {
    degree = 89;
  } else if(field == Field::GF2_79) {
    degree = 79;
  } else {
    throw std::invalid_argument("Unsupported binary field");
  }

  uint131_t result = make_uint131(1);
  for(int i = degree - 1; i >= 0; i--) {
    result = mul(result, result, field);
    if(i != 0)
      result = mul(result, a, field);
  }
  return result;
}

} // namespace gf2

#endif
