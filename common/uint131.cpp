#include <cassert>
#include <random>
#include <stdio.h>
#include <string>
#include <vector>

#include "uint131.h"

namespace {

uint8_t from_hex(char hex)
{
  if(hex >= 'a' && hex <= 'f') {
      return hex - 'a' + 10;
  } else if(hex >= 'A' && hex <= 'F') {
      return hex - 'A' + 10;
  } else if(hex >= '0' && hex <= '9') {
      return hex - '0';
  }

  return '0';
}

std::vector<uint64_t> parse_hex_uint64(const std::string& s)
{
    std::string hex = s;

    // Pad with leading zeros
    if (hex.length() % 16 != 0)
    {
      hex = std::string(16 - hex.length() % 16, '0') + hex;
    }

    std::vector<uint64_t> ara;

    for (int i = 0; i < hex.length() / 16; i++)
    {
      uint64_t word = 0;
      for (int j = 0; j < 8; j++)
      {
        word <<= 8;
        uint8_t high = from_hex(hex[16 * i + 2 * j]);
        uint8_t low = from_hex(hex[16 * i + 2 * j + 1]);

        uint8_t value = (high << 4) | low;
        word |= value;
      }
      ara.push_back(word);
    }

    return ara;
}

}

uint131_t make_uint131(uint32_t x)
{
  uint131_t val;
  memset(&val, 0, sizeof(val));

  val.v[0] = x;

  return val;
}

uint131_t make_uint131(const std::string& hex)
{

  assert(hex.length() <= 33);

  uint131_t val;
  memset(&val, 0, sizeof(val));

  std::vector<uint64_t> words = parse_hex_uint64(hex);

  for(int i = 0; i < words.size(); i++) {
    val.v[i] = words[words.size() - 1 - i];
  }
  
  return val;
}

bool operator==(const uint131_t& a, const uint131_t& b)
{
  for(int i = 0; i < 3; i++) {
    if(a.v[i] != b.v[i]) {
      return false;
    }
  }

  return true;
}

bool operator!=(const uint131_t& a, const uint131_t& b)
{
  return !(a == b);
}

bool is_odd(const uint131_t& x)
{
  return x.v[0] & 0x01;
}

std::string to_str(const uint131_t& x)
{
  char buf[256] = "";

  sprintf(buf, "%.1lX%.16lX%.16lX", x.v[2], x.v[1], x.v[0]);

  return std::string(buf);
}
