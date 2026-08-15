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

inline void square(const uint131_t& a, uint64_t product[6])
{
  const uint64_t words[3] = {a.w.v0, a.w.v1, a.w.v2};
  for(int i = 0; i < 3; i++) {
    clmul64(words[i], words[i], product[2 * i], product[2 * i + 1]);
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

inline uint131_t square(const uint131_t& a, Field field)
{
  uint64_t product[6] = {};
  detail::square(a, product);

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

namespace detail {

inline uint131_t square_n(uint131_t value, int count, Field field)
{
  for(int i = 0; i < count; i++) value = square(value, field);
  return value;
}

inline uint131_t inv_131(const uint131_t& x)
{
  // Optimal addition chain for 130: 1,2,4,8,16,32,64,128,130.
  uint131_t x2 = mul(square_n(x, 1, Field::GF2_131), x, Field::GF2_131);
  uint131_t x4 = mul(square_n(x2, 2, Field::GF2_131), x2, Field::GF2_131);
  uint131_t x8 = mul(square_n(x4, 4, Field::GF2_131), x4, Field::GF2_131);
  uint131_t x16 = mul(square_n(x8, 8, Field::GF2_131), x8, Field::GF2_131);
  uint131_t x32 = mul(square_n(x16, 16, Field::GF2_131), x16, Field::GF2_131);
  uint131_t x64 = mul(square_n(x32, 32, Field::GF2_131), x32, Field::GF2_131);
  uint131_t x128 = mul(square_n(x64, 64, Field::GF2_131), x64, Field::GF2_131);
  uint131_t x130 = mul(square_n(x128, 2, Field::GF2_131), x2, Field::GF2_131);
  return square(x130, Field::GF2_131);
}

inline uint131_t inv_89(const uint131_t& x)
{
  // Optimal addition chain for 88: 1,2,4,8,16,32,64,80,88.
  uint131_t x2 = mul(square_n(x, 1, Field::GF2_89), x, Field::GF2_89);
  uint131_t x4 = mul(square_n(x2, 2, Field::GF2_89), x2, Field::GF2_89);
  uint131_t x8 = mul(square_n(x4, 4, Field::GF2_89), x4, Field::GF2_89);
  uint131_t x16 = mul(square_n(x8, 8, Field::GF2_89), x8, Field::GF2_89);
  uint131_t x32 = mul(square_n(x16, 16, Field::GF2_89), x16, Field::GF2_89);
  uint131_t x64 = mul(square_n(x32, 32, Field::GF2_89), x32, Field::GF2_89);
  uint131_t x80 = mul(square_n(x64, 16, Field::GF2_89), x16, Field::GF2_89);
  uint131_t x88 = mul(square_n(x80, 8, Field::GF2_89), x8, Field::GF2_89);
  return square(x88, Field::GF2_89);
}

inline uint131_t inv_79(const uint131_t& x)
{
  // Optimal addition chain for 78: 1,2,3,6,12,24,48,72,78.
  uint131_t x2 = mul(square_n(x, 1, Field::GF2_79), x, Field::GF2_79);
  uint131_t x3 = mul(square_n(x2, 1, Field::GF2_79), x, Field::GF2_79);
  uint131_t x6 = mul(square_n(x3, 3, Field::GF2_79), x3, Field::GF2_79);
  uint131_t x12 = mul(square_n(x6, 6, Field::GF2_79), x6, Field::GF2_79);
  uint131_t x24 = mul(square_n(x12, 12, Field::GF2_79), x12, Field::GF2_79);
  uint131_t x48 = mul(square_n(x24, 24, Field::GF2_79), x24, Field::GF2_79);
  uint131_t x72 = mul(square_n(x48, 24, Field::GF2_79), x24, Field::GF2_79);
  uint131_t x78 = mul(square_n(x72, 6, Field::GF2_79), x6, Field::GF2_79);
  return square(x78, Field::GF2_79);
}

} // namespace detail

inline uint131_t inv(const uint131_t& a, Field field)
{
  if(a == make_uint131(0)) return {};

  switch(field) {
  case Field::GF2_131:
    return detail::inv_131(a);
  case Field::GF2_89:
    return detail::inv_89(a);
  case Field::GF2_79:
    return detail::inv_79(a);
  default:
    throw std::invalid_argument("Unsupported binary field");
  }
}

} // namespace gf2

#endif
