#ifndef _EC_RHO_H
#define _EC_RHO_H

#include <functional>
#include <iostream>

#include "ecc.h"
#include "fmt/format.h"
#include "montgomery.h"

struct DistinguishedPoint {
  uint131_t a;
  ecc::ecpoint_t p;
  uint64_t length = 0;
  int dp_bits = 0;

  DistinguishedPoint() {}

  DistinguishedPoint(uint131_t a, ecc::ecpoint_t p, int dp_bits, uint64_t length)
      : a(a), p(p), dp_bits(dp_bits)
  {
    this->length = length;
  }
};

struct DPData {
  uint8_t sign;
  uint8_t a[17];
};

struct EncodedDP {
  // This will be the key in the database
  uint8_t tx[17];

  // This will be the value in the database
  DPData data;

  // Walk length
  uint8_t len[5];

  // A few bits of the Y coordinate.
  uint8_t checksum;

  EncodedDP() { memset(tx, 0, sizeof(tx)); }
};

struct DPHeader {
  uint8_t version;
  uint32_t count;
  uint8_t curve_id;
  int dp_bits;
};

struct RWPoint {
  ecc::ecpoint_t p;
  uint131_t a;
  uint131_t b;
};

class DistinguishedPointFinder {

public:
  virtual size_t work_per_step() = 0;
  virtual int iters_per_step() = 0;
  virtual void
  set_callback(std::function<void(const std::vector<DistinguishedPoint>& p)> callback) = 0;
  virtual void init() = 0;
  virtual void init(const std::string& file) = 0;
  virtual double step() = 0;
  virtual void save_progress(const std::string& file) = 0;
  virtual int parallel_walks() = 0;
  virtual ~DistinguishedPointFinder() {}
};

std::vector<RWPoint> get_rw_points();

std::vector<DistinguishedPoint> decode_dps(const uint8_t* bytes, size_t size, bool verify = false);
DistinguishedPoint decode_dp(const EncodedDP& dp, int dpbits, bool verify = false);
std::vector<uint8_t> encode_dps(const std::vector<DistinguishedPoint>& dps, int curve, int dpbits);
bool verify_dp(const DistinguishedPoint& dp);

class RhoSolver {

private:
  std::vector<RWPoint> _rwpoints = get_rw_points();
  ecc::ecpoint_t _dp;
  uint131_t _a1;
  uint131_t _a2;

public:
  RhoSolver(const ecc::ecpoint_t dp, uint131_t a1, uint131_t a2) : _dp(dp), _a1(a1), _a2(a2) {}

  RWPoint iterate(RWPoint p)
  {
    int idx = p.p.x.v[0] & 0x1f;

    p.a = mont::add_mod_n(p.a, _rwpoints[idx].a, ecc::n());
    p.b = mont::add_mod_n(p.b, _rwpoints[idx].b, ecc::n());

    p.p = ecc::add(p.p, _rwpoints[idx].p);

    return p;
  }

  uint64_t get_length(uint131_t a, uint64_t max_len)
  {
    RWPoint p;
    p.a = a;
    p.p = ecc::mul(a, ecc::g());
    uint131_t b = make_uint131(0);

    uint64_t count = 0;

    while(p.p.x != _dp.x && count < max_len) {
      p = iterate(p);
      count++;
    }

    return count;
  }

  void solve()
  {
    std::cout << "Getting walk 1 length" << std::endl;
    uint64_t len1 = get_length(_a1, 0xffffffff);
    std::cout << len1 << std::endl;
    uint64_t len2 = get_length(_a2, 0xffffffff);
    std::cout << len2 << std::endl;

    if(len1 < len2) {
      std::swap(len1, len2);
      std::swap(_a1, _a2);
    }

    // Starting points only have the 'a' exponent. b = 0
    ecc::ecpoint_t p1 = ecc::mul(_a1, ecc::g());
    uint131_t b1 = make_uint131(0);

    ecc::ecpoint_t p2 = ecc::mul(_a2, ecc::g());
    uint131_t b2 = make_uint131(0);

    std::cout << fmt::format("x={0}    y={1}", to_str(p1.x), to_str(p2.x)) << std::endl;
    std::cout << fmt::format("x={0}    y={1}", to_str(p1.y), to_str(p2.y)) << std::endl;

    uint64_t diff = len1 - len2;

    // Advance the longer walk until both walks have same remaining steps
    RWPoint rw1;
    rw1.a = _a1;
    rw1.b = make_uint131(0);
    rw1.p = p1;

    for(uint64_t i = 0; i < diff; i++) {
      rw1 = iterate(rw1);
    }

    RWPoint rw2;
    rw2.a = _a2;
    rw2.b = make_uint131(0);
    rw2.p = p2;

    std::cout << "Starting points:" << std::endl;
    std::cout << fmt::format("x={0}     x={1}", to_str(p1.x), to_str(p2.x)) << std::endl;
    std::cout << fmt::format("x={0}     x={1}", to_str(p1.y), to_str(p2.y)) << std::endl;

    // Advance both walks until collision
    uint64_t remaining = len2;
    while(rw1.p.x != rw2.p.x && remaining > 0) {
      rw1 = iterate(rw1);
      rw2 = iterate(rw2);
      remaining--;
    }

    std::cout << "Colliding point:" << std::endl;
    std::cout << fmt::format("x={0}   x={1}", to_str(rw1.p.x), to_str(rw2.p.x)) << std::endl;
    std::cout << fmt::format("y={0}   y={1}", to_str(rw1.p.y), to_str(rw2.p.y)) << std::endl;
    std::cout << fmt::format("a={0}   a={1}", to_str(rw1.a), to_str(rw2.a)) << std::endl;
    std::cout << fmt::format("b={0}   b={1}", to_str(rw1.b), to_str(rw2.b)) << std::endl;

    // Calculate the discrete logarithm
    // a1G + b1Q = a2G + b2Q
    // (a1 - a2)G = (b2 - b1)Q
    // (a1 - a2)G = (b2 - b1)kG
    // (a1 - a2) = (b2 - b1)k
    // (a1 - a2) / (b2 - b1) = k
    uint131_t d1 = mont::sub_mod_n(rw1.a, rw2.a, ecc::n());
    uint131_t d2 = mont::sub_mod_n(rw2.b, rw1.b, ecc::n());
    uint131_t k = mont::mul_mod_n(d1, mont::inv_mod_n(d2, ecc::n()), ecc::n());

    // Validate
    ecc::ecpoint_t result_q = ecc::mul(k, ecc::g());
    ecc::ecpoint_t real_q = ecc::q();

    if(!ecc::is_equal(result_q, real_q)) {
      std::cout << "Something went wrong" << std::endl;
      std::cout << fmt::format("Expected {} {}", to_str(real_q.x), to_str(real_q.y)) << std::endl;
      std::cout << fmt::format("Got      {} {}", to_str(result_q.x), to_str(result_q.y))
                << std::endl;
    } else {
      std::cout << fmt::format("k={}", to_str(k)) << std::endl;
    }
  }
};

#endif