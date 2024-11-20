#ifndef _UINT131_H
#define _UINT131_H

#include <stdint.h>
#include <string.h>
#include <string>

struct uint131_t {
  uint64_t v[3];
};

uint131_t make_uint131(uint32_t x);
uint131_t make_uint131(const std::string& hex);

bool operator==(const uint131_t& a, const uint131_t& b);
bool operator!=(const uint131_t& a, const uint131_t& b);
bool is_odd(const uint131_t& x);

std::string to_str(const uint131_t& x);

#endif