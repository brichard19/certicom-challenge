#ifndef _EC_RHO_H
#define _EC_RHO_H

#include <functional>

#include "ecc.h"

struct DistinguishedPoint {
  uint131_t a;
  ecc::ecpoint_t p;
  uint64_t length = 0;
  int dp_bits = 0;

  DistinguishedPoint(){}

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

#endif