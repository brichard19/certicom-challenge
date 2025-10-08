#include <random>

#include "ecc.h"
#include "ecc_internal.h"
#include "montgomery.h"
#include "util.h"
#include <stdexcept>


CurveParameters _eccp131 = {
    .p = {{0x194c43186b3abc0b, 0x8e1d43f293469e33, 0x4}},
    .a = {{0x9893f47bd114555d, 0x2b935e3727d799fa, 0x1}},
    .b = {{0xe59d7a40bc8648b, 0x55abb9bc04daa824, 0x4}},
    .n = {{0x7f7ed728f6b8e6f1, 0x8e1d43f293469e31, 0x4}},

    .gx = {{0x11106cddc12294ab, 0x862b445b2d90f966, 0x1}},
    .gy = {{0xbdd264b969654eb3, 0x2ccf5a5df57b1433, 0x3}},

    .qx = {{0xa036d1676c1cec9c, 0x9bf4f50b596c4922, 0x1}},

    .qy = {{0xc9a258ee3a5aa43f, 0xfd3113f721b1e7aa, 0x3}},

    // k such that k * k^-1 = -1 (mod R)
    // k = (r * r_inv -1) // p
    .k = {{0xe0587d72985b105d, 0xf1fd54b0309e1ab9, 0x7}},

    .one = {{0xe6b3bce794c543f5, 0x71e2bc0d6cb961cc, 0x3}},

    .p_minus_2 = {{0x194c43186b3abc09, 0x8e1d43f293469e33, 0x4}},

    // (P + 1) / 4
    .sqrt = {{0xc65310c61aceaf03, 0x238750fca4d1a78c, 0x1}},

    // R mod p
    .r = {{0xe6b3bce794c543f5, 0x71e2bc0d6cb961cc, 0x3}},

    // R^2
    .r2 = {{0xa0a3d1806d7fee3d, 0x93a974dd09a9f553, 0x1}},

    .bits = 131,

    .words = (131 + 63) / 64,

    .name = "ecp131"
};

CurveParameters _eccp79 = {
  .p = {{0x5177412aca899cf5, 0x62ce, 0x00}},

  .a = {{0x1274945cc02f5789, 0x5320, 0x0}},

  .b = {{0x967af723a2d4326f, 0x1e5d, 0x0}},

  .n = {{0x5177407b7258dc31, 0x62ce, 0x00}},

  .gx = {{0x5e14178a38bb98f0, 0x3088, 0x0}},

  .gy = {{0x9612fecbf0419185, 0x4f96, 0x0}},

  .qx = {{0x4ee59e29999d86f, 0x4e23, 0x0}},

  .qy = {{0xc994049f2ed05e7f, 0x2f1d, 0x0}},

  .k = {{0x6e3655426732d0a3, 0xcafea9fd045a89b6, 0x4}},

  .one = {{0x150619962eb12300, 0x5b4b, 0x0}},

  .p_minus_2 = {{0x5177412aca899cf3, 0x62ce, 0x00}},

  // (p - 5) / 8
  .sqrt = {{0xca2ee8255951339e, 0xc59, 0x00}},
  
  .r = {{0x150619962eb12300, 0x5b4b, 0x0}},

  // R^2
  .r2 = {{0x45cc4bd33e72e43c, 0x1b44, 0x0}},

  .bits = 79,

  .words = (79 + 63) / 64,

  .name = "ecp79"
};

CurveParameters _params;

namespace ecc {
void set_curve(const std::string& curve_name)
{
  if(curve_name == "ecp131") {
    _params = _eccp131;
  } else if(curve_name == "ecp79") {
    _params = _eccp79;
  } else {
    throw std::runtime_error("Invalid curve name");
  }
}

}

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
  uint131_t ax = mont::mul(_params.a, p.x);

  uint131_t rs = mont::add(mont::add(x3, ax), _params.b);

  return y2 == rs;
}

ecpoint_t dbl(const ecpoint_t& p)
{
  if(is_infinity(p)) {
    return p;
  }

  // 3x^2 + a / 2y

  uint131_t x2 = mont::square(p.x);
  uint131_t rise = mont::add(mont::add(mont::add(x2, x2), x2), _params.a);

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
    if(_params.p.w.v2 == 0) {
      r.w.v2 = 0;
      r.w.v1 %= (_params.p.w.v1 + 1);

      gte = r.w.v0 >= _params.p.w.v0;
    } else {
      r.w.v2 %= (_params.p.w.v2 + 1);
      gte = r.w.v1 >= _params.p.w.v1;
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
  std::vector<ecpoint_t> qmul(_params.bits + 1);

  // P = kQ
  std::vector<ecpoint_t> p(k.size());

  qmul[0] = q;

  for(int i = 1; i < qmul.size(); i++) {
    qmul[i] = ecc::dbl(qmul[i - 1]);
  }

  // For each bit in the private keys
  for(int b = 0; b < _params.bits; b++) {

    // Buffer to hold the multiplication chain
    std::vector<uint131_t> mbuf(k.size());

    uint131_t t;
    uint131_t inverse = _params.one;

    for(int i = 0; i < k.size(); i++) {

      int bit = get_bit(k[i], b);

      if(bit == 0 || is_infinity(p[i]) || is_equal(p[i], qmul[b])) {
        // Nothing to do. Multiply by 1
        t = _params.one;
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
  return ecpoint_t(_params.gx, _params.gy);
}

ecpoint_t q()
{
  return ecpoint_t(_params.qx, _params.qy);
}

uint131_t p()
{
  return _params.p;
}

uint131_t a()
{
  return _params.a;
}

uint131_t b()
{
  return _params.b;
}

uint131_t n()
{
  return _params.n;
}

std::string curve_name()
{
  return _params.name;
}

int curve_strength()
{
  return _params.bits;
}

}