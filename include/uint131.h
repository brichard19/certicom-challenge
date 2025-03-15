#ifndef _UINT131_H
#define _UINT131_H

#include <stdint.h>
#include <string.h>
#include <string>

union uint131_t {
  uint64_t v[3];
  struct {
    uint64_t v0;
    uint64_t v1;
    uint32_t v2;
    uint32_t padding;
  }w;
};

union uint160_t {
  uint64_t v[3];
  struct {
    uint64_t v0;
    uint64_t v1;
    uint32_t v2;
    uint32_t padding;
  }w;
};

struct uint262_t {
  uint64_t v[5];
};

uint131_t make_uint131(uint32_t x);
uint131_t make_uint131(const std::string& hex);

bool operator==(const uint131_t& a, const uint131_t& b);
bool operator!=(const uint131_t& a, const uint131_t& b);
bool is_odd(const uint131_t& x);

std::string to_str(const uint131_t& x);

#endif