#ifndef _EC_RHO_H
#define _EC_RHO_H

#include <functional>

#include "ecc.h"

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

};


std::vector<RWPoint> get_rw_points();

std::vector<DistinguishedPoint> decode_dps(const uint8_t* bytes, size_t size);
std::vector<uint8_t> encode_dps(const std::vector<DistinguishedPoint>& dps, int dbits, int curve);

#endif