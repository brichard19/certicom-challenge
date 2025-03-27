#ifndef _ECC_H
#define _ECC_H

#include <vector>

#include "uint131.h"

namespace ecc {

void set_curve(const std::string& curve_name);

struct ecpoint_t {

  uint131_t x;
  uint131_t y;

  ecpoint_t()
  {
    this->x = make_uint131((uint32_t)0);
    this->x.w.v2 = (uint32_t)-1;
    this->y = make_uint131((uint32_t)0);
  }

  ecpoint_t(const uint131_t& x, const uint131_t& y)
  {
    this->x = x;
    this->y = y;
  }
};

bool is_infinity(const ecpoint_t& p);
bool is_equal(const ecpoint_t& p, const ecpoint_t& q);
bool is_neg(const ecpoint_t& p, const ecpoint_t& q);
bool exists(const ecpoint_t& p);

ecpoint_t dbl(const ecpoint_t& p);
ecpoint_t add(const ecpoint_t& p, const ecpoint_t& q);
ecpoint_t mul(const uint131_t& k, const ecpoint_t& p);
ecpoint_t g();
ecpoint_t q();
uint131_t p();
uint131_t a();
uint131_t b();
uint131_t n();
std::string curve_name();
int curve_strength();

uint131_t genkey();
std::vector<uint131_t> genkeys(int count, int seed);

}

#endif