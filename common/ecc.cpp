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