#ifndef _EC_RHO_H
#define _EC_RHO_H

#include <functional>

#include "ecc.h"

#define DP_BITS 25
#define X_TRUNC_LEN (((131 - DP_BITS) + 7) / 8)

// TX = Truncated x coordinate (no distinguished bits)
// S = sign bit (1 byte)
// A = a exponent
// L = walk length (5 bytes)
//
// [        TX        ][S][        A        ][    L    ]
// Encodes a DistinguishedPoint into a string of bytes

#define X_OFFSET 0
#define SIGN_OFFSET X_TRUNC_LEN
#define A_OFFSET (SIGN_OFFSET + 1)
#define LEN_OFFSET (A_OFFSET + 17)
#define ENCODED_DP_SIZE (LEN_OFFSET + 5)

struct DistinguishedPoint {
  uint131_t a;
  ecc::ecpoint_t p;
  uint64_t length = 0;

  DistinguishedPoint()
  {
  }

  DistinguishedPoint(uint131_t a, ecc::ecpoint_t p, uint64_t length=0) : a(a), p(p)
  {
    this->length = length;
  }
};

struct DPHeader {
  uint8_t version;
  uint32_t count;
  uint8_t curve;
  uint8_t dbits;
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
  virtual void set_callback(std::function<void(const std::vector<DistinguishedPoint>& p)> callback) = 0;
  virtual void init() = 0;
  virtual void init(const std::string& file) = 0;
  virtual double step() = 0;
  virtual void save_progress(const std::string& file) = 0;
  virtual int parallel_walks() = 0;
  virtual ~DistinguishedPointFinder(){}
};


std::vector<RWPoint> get_rw_points();

std::vector<DistinguishedPoint> decode_dps(const uint8_t* bytes, size_t size);
DistinguishedPoint decode_dp(const uint8_t* bytes, int dbits);
std::vector<uint8_t> encode_dps(const std::vector<DistinguishedPoint>& dps, int dbits, int curve);

#endif