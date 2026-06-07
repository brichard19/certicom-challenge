#ifndef _RNG_h
#define _RNG_H

#include <random>
#include <stdint.h>

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
  IntRNG() { _gen = std::mt19937(_rd()); }

  virtual uint64_t next() { return _d(_gen); }
};

class DeterministicRNG : public RNG {

private:
  uint64_t _state;

public:
  DeterministicRNG(int seed) { _state = seed; }

  virtual uint64_t next()
  {
    _state = _state * 6364136223846793005 + 1442695040888963407;

    return _state;
  }
};

#endif