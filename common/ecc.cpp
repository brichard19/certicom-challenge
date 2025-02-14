#include <random>

#include "ecc.h"
#include "eccp131.h"
#include "mont.h"
#include "util.h"

namespace {

class IntRNG {

private:
  std::random_device _rd;
  std::mt19937 _gen;
  std::uniform_int_distribution<uint64_t> _d;

public:
  IntRNG()
  {
    _gen = std::mt19937(_rd());
  }

  uint64_t next()
  {
    return _d(_gen);
  }
};

// Internal RNG
IntRNG _rng;


class DeterministicRNG{

private:
    uint64_t _state;

public:
    DeterministicRNG(int seed)
    {
        _state = seed;
    }

    uint64_t next()
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

  p.x.v[2] = (uint64_t)-1;

  return p;
}

bool is_infinity(const ecpoint_t& p)
{
  // If the high bits of x are all 1, then its a point-at-infinity
  return p.x.v[2] == (uint64_t)-1;
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

  if(p.y == mont::neg(q.y)) {
    return false;
  }

  return true;
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

uint131_t genkey()
{
  uint131_t r;

  bool gte = true;

  memset(&r, 0, sizeof(r));

  do {
    for(int i = 0; i < _words; i++) {
      r.v[i] = _rng.next();
    }

    // Limit most-significant word
    r.v[_words - 1] = r.v[_words - 1] % (_n.v[_words - 1] + 1);

    // Check if less than P
    for(int i = _words - 1; i >= 0; i--) {
      if(r.v[i] < _n.v[i]) {
        gte = false;
        break;
      } else if(r.v[i] > _n.v[i]) {
        gte = true;
        break;
      }
    }
  }while (gte);

  return r;
}

// Use this when you need to generate keys in a deterministic way i.e. the same seed will
// always generate the same keys
std::vector<uint131_t> genkeys(int count, int seed)
{
  DeterministicRNG rng(seed);

  std::vector<uint131_t> keys;

  for(int k = 0; k < count; k++) {
    uint131_t r;
    memset(&r, 0, sizeof(r));

    bool gte = true;

    do {
      for(int i = 0; i < _words; i++) {
        r.v[i] = rng.next();
      }

      // Limit most-significant word
      r.v[_words - 1] = r.v[_words - 1] % (_n.v[_words - 1] + 1);

      // Check if less than P
      for(int i = _words - 1; i >= 0; i--) {
        if(r.v[i] < _n.v[i]) {
          gte = false;
          break;
        } else if(r.v[i] > _n.v[i]) {
          gte = true;
          break;
        }
      }
    }while (gte);

    keys.push_back(r);
  }

  return keys;
}

ecpoint_t mul(const uint131_t& k, const ecpoint_t& p)
{

  ecc::ecpoint_t sum;

  ecc::ecpoint_t adder = p;

  uint64_t mask = 0x01;
  for(int i = 0; i < 2; i++) {

    mask = 0x01;
    for(int bit = 0; bit < 64; bit++) {
      if(k.v[i] & mask) {
        sum = ecc::add(sum, adder);
      }

      mask <<= 1;
      adder = ecc::dbl(adder);
    }
  }

  mask = 0x01;
  for(int bit = 0; bit < 3; bit++) {
    if(k.v[2] & mask ) {
      sum = ecc::add(sum, adder);
    }

    mask <<= 1;
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