#include <random>

#include "ecc.h"
#include "eccp131.h"
#include "mont.h"
#include "util.h"

namespace {

class RNG {

public:
  virtual uint64_t next() = 0;
};

class IntRNG : public RNG {

private:
  std::random_device _rd;
  std::mt19937 _gen;
  std::uniform_int_distribution<uint64_t> _d;

public:
  IntRNG()
  {
    _gen = std::mt19937(_rd());
  }

  virtual uint64_t next()
  {
    return _d(_gen);
  }
};

// Internal RNG
IntRNG _rng;


class DeterministicRNG : public RNG {

private:
    uint64_t _state;

public:
    DeterministicRNG(int seed)
    {
        _state = seed;
    }

    virtual uint64_t next()
    {
        _state = _state * 6364136223846793005 + 1442695040888963407;

        return _state;
    }
};

}

namespace ecc {


ecpoint_t make_infinity()
{
  ecpoint_t p;

  memset(&p, 0, sizeof(p));

  p.x.w.v2 = (uint32_t)-1;
  return p;
}

bool is_infinity(const ecpoint_t& p)
{
  // If the high bits of x are all 1, then its a point-at-infinity
  return p.x.w.v2 == (uint32_t)-1;
}

bool is_equal(const ecpoint_t& p, const ecpoint_t& q)
{
  return p.x == q.x && p.y == q.y;
}

bool is_neg(const ecpoint_t& p, const ecpoint_t& q)
{
  // P == -Q when Px == Qx and Py == -Qy
  if(p.x != q.x) {
    return false;
  }

  return p.y == mont::neg(q.y);
}

bool exists(const ecpoint_t& p)
{
  uint131_t y2 = mont::square(p.y);
  uint131_t x3 = mont::mul(p.x, mont::square(p.x));
  uint131_t ax = mont::mul(_a, p.x);

  uint131_t rs = mont::add(mont::add(x3, ax), _b);

  return y2 == rs;
}

ecpoint_t dbl(const ecpoint_t& p)
{
  if(is_infinity(p)) {
    return p;
  }

  // 3x^2 + a / 2y

  uint131_t x2 = mont::square(p.x);
  uint131_t rise = mont::add(mont::add(mont::add(x2, x2), x2), _a);

  uint131_t run = mont::add(p.y, p.y);
  uint131_t s = mont::mul(rise, mont::inv(run));

  uint131_t s2 = mont::square(s);

  uint131_t tmp1 = mont::sub(s2, p.x);
  uint131_t x = mont::sub(tmp1, p.x);

  uint131_t tmp2 = mont::sub(p.x, x);
  uint131_t tmp3 = mont::mul(s, tmp2);
  uint131_t y = mont::sub(tmp3, p.y);

  return ecc::ecpoint_t(x, y);

}

ecpoint_t add(const ecpoint_t& p, const ecpoint_t& q)
{
  if(is_infinity(p)) {
    return q;
  }

  if(is_infinity(q)) {
    return p;
  }
  
  if(is_equal(p, q)) {
    return dbl(p);
  }
  
  if(is_neg(p, q)) {
    return make_infinity();
  }

  uint131_t rise = mont::sub(p.y, q.y);
  uint131_t run = mont::sub(p.x, q.x);

  // s = inv(rise / run)
  uint131_t s = mont::inv(run);
  s = mont::mul(s, rise);

  uint131_t s2 = mont::square(s);

  uint131_t tmp1 = mont::sub(s2, p.x);
  uint131_t x = mont::sub(tmp1, q.x);

  uint131_t tmp2 = mont::sub(p.x, x);
  uint131_t tmp3 = mont::mul(s, tmp2);
  uint131_t y = mont::sub(tmp3, p.y);

  return ecc::ecpoint_t(x, y);
}

uint131_t genkey(RNG& rng)
{
  uint131_t r;

  bool gte = true;

  memset(&r, 0, sizeof(r));

  do {
    r.w.v0 = rng.next();
    r.w.v1 = rng.next();
    r.w.v2 = (uint32_t)rng.next();

    // mod p
    if(_p.w.v2 == 0) {
      r.w.v2 = 0;
      r.w.v1 %= (_p.w.v1 + 1);

      gte = r.w.v0 >= _p.w.v0;
    } else {
      r.w.v2 %= (_p.w.v2 + 1);
      gte = r.w.v1 >= _p.w.v1;
    }

  }while (gte);

  return r;
}

uint131_t genkey()
{
  return genkey(_rng);
}

// Use this when you need to generate keys in a deterministic way i.e. the same seed will
// always generate the same keys
std::vector<uint131_t> genkeys(int count, int seed)
{
  DeterministicRNG rng(seed);

  std::vector<uint131_t> keys;

  for(int k = 0; k < count; k++) {
    uint131_t r;
    r = genkey(rng);

    keys.push_back(r);
  }

  return keys;
}

ecpoint_t mul(const uint131_t& k, const ecpoint_t& p)
{
  ecc::ecpoint_t sum;

  ecc::ecpoint_t adder = p;

  // v0
  uint64_t mask = 0x01;
  for(int i = 0; i < 64; i++) {
      if(k.w.v0 & mask) {
        sum = ecc::add(sum, adder);
      }

      mask <<= 1;
      adder = ecc::dbl(adder);
  }

  // v1
  mask = 0x01;
  for(int i = 0; i < 64; i++) {
      if(k.w.v1 & mask) {
        sum = ecc::add(sum, adder);
      }

      mask <<= 1;
      adder = ecc::dbl(adder);
  } 

  // v2
  uint32_t mask32 = 0x01;
  for(int i = 0; i < 32; i++) {
      if(k.w.v2 & mask32) {
        sum = ecc::add(sum, adder);
      }

      mask32 <<= 1;
      adder = ecc::dbl(adder);
  } 

  return sum;
}

int get_bit(uint131_t x, int bit)
{
  if(bit >= 128) {
    return (x.w.v2 >> (bit & 0x7f)) & 1;
  } else if(bit >= 64) {
    return (x.w.v1 >> (bit & 0x3f)) & 1;
  } else {
    return (x.w.v0 >> bit) & 1;
  }
}

std::vector<ecpoint_t> mul(const std::vector<uint131_t>& k, const ecpoint_t& q)
{
  // Create lookup table for Q, 2Q, 4Q... 2^131Q
  std::vector<ecpoint_t> qmul(_bits + 1);

  // P = kQ
  std::vector<ecpoint_t> p(k.size());

  qmul[0] = q;

  for(int i = 1; i < qmul.size(); i++) {
    qmul[i] = ecc::dbl(qmul[i - 1]);
  }

  // For each bit in the private keys
  for(int b = 0; b < _bits; b++) {

    // Buffer to hold the multiplication chain
    std::vector<uint131_t> mbuf(k.size());

    uint131_t t;
    uint131_t inverse = _one;

    for(int i = 0; i < k.size(); i++) {

      int bit = get_bit(k[i], b);

      if(bit == 0 || is_infinity(p[i]) || is_equal(p[i], qmul[b])) {
        // Nothing to do. Multiply by 1
        t = _one;
      } else {
        // x2 - x1
        t = mont::sub(qmul[b].x, p[i].x);
      }
      inverse = mont::mul(inverse, t);
      mbuf[i] = inverse;
    }

    inverse = mont::inv(inverse);

    for(int i = k.size() - 1; i >= 0; i--) {

      uint131_t px = p[i].x;
      uint131_t py = p[i].y;

      int bit = get_bit(k[i], b);

      if(!bit) {
        // Nothing to do
        continue;
      } else if(is_infinity(p[i])) {
        p[i] = qmul[b];
        continue;
      } else if(is_equal(p[i], qmul[b])) {
        // Double
        p[i] = qmul[b + 1];
        continue;
      }

      uint131_t s;
      if(i > 0) {
        // Get the 2nd-last element (product of all factors up to that number)
        // e.g. abcd
        uint131_t m = mbuf[i - 1];

        // Multiply to cancel out all factors except the last one
        // e.g. abcd * (abcde)^-1 = e^-1
        s = mont::mul(inverse, m);

        // Cancel out from the inverse
        // e.g. abcde * e^-1 = abcd
        uint131_t diff = mont::sub(qmul[b].x, p[i].x);
        
        inverse = mont::mul(inverse, diff);
      } else {
        s = inverse;
      }
      // Complete addition

      const uint131_t qx = qmul[b].x;
      const uint131_t qy = qmul[b].y;

      uint131_t rise = mont::sub(qy, py);
      s = mont::mul(s, rise);
      uint131_t s2 = mont::square(s);

      uint131_t tmp1 = mont::sub(s2, px);
      uint131_t x = mont::sub(tmp1, qx);

      uint131_t tmp2 = mont::sub(px, x);
      uint131_t tmp3 = mont::mul(s, tmp2);
      uint131_t y = mont::sub(tmp3, py);

      p[i].x = x;
      p[i].y = y;

    }
  }

  return p;
}

// Calculate y from x and sign
uint131_t calc_y(const uint131_t& x, int sign)
{
  // y^2 = x^3 + ax + b
  // y = sqrt(x^3 + ax + b)

  uint131_t x3 = mont::mul(mont::square(x), x);
  uint131_t ax = mont::mul(ecc::a(), x);

  uint131_t y2 = mont::add(mont::add(x3, ax), ecc::b());

  uint131_t y = mont::sqrt(y2);

  if(get_bit(y, 0) == sign) {
    return y;
  }
  
  return mont::sub(ecc::p(), y);
}

ecpoint_t g()
{
  return ecpoint_t(_gx, _gy);
}

ecpoint_t q()
{
  return ecpoint_t(_qx, _qy);
}

uint131_t p()
{
  return _p;
}

uint131_t a()
{
  return _a;
}

uint131_t b()
{
  return _b;
}

uint131_t n()
{
  return _n;
}

std::string curve_name()
{
  return _curve_name;
}

int curve_strength()
{
  return _bits;
}

}