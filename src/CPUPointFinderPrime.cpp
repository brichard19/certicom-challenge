#include "CPUPointFinderPrime.h"

#include <cassert>
#include <chrono>
#include <vector>

#include "ecc_internal.h"
#include "montgomery.h"

CPUPointFinderPrime::CPUPointFinderPrime(int dpbits, size_t num_points)
    : CPUPointFinder(dpbits, num_points)
{
}

void CPUPointFinderPrime::do_step_prime(std::vector<DistinguishedPoint>& results)
{
  for(size_t i = 0; i < _num_points; i++) {
    if((_x[i].w.v0 & _dpmask) == 0) {
      ecc::ecpoint_t point(_x[i], _y[i]);
      assert(ecc::exists(point));
      results.emplace_back(_priv[i], point, _dpbits, _walk_len[i]);
      generate_walk(i);
    }

    size_t rw_index = _x[i].w.v0 & 0x1f;
    uint131_t denominator = mont::sub(_rx[rw_index], _x[i]);
    _chain[i] = mont::mul(i == 0 ? _params.one : _chain[i - 1], denominator);
  }

  uint131_t inverse = mont::inv(_chain.empty() ? _params.one : _chain.back());

  for(size_t i = _num_points; i-- > 0;) {
    size_t rw_index = _x[i].w.v0 & 0x1f;
    ecc::ecpoint_t p(_x[i], _y[i]);
    ecc::ecpoint_t q(_rx[rw_index], _ry[rw_index]);

    uint131_t denominator = mont::sub(q.x, p.x);
    uint131_t denominator_inverse = i == 0 ? inverse : mont::mul(inverse, i >= 1 ? _chain[i - 1] : _params.one);
    inverse = mont::mul(inverse, denominator);

    uint131_t slope = mont::mul(mont::sub(q.y, p.y), denominator_inverse);
    uint131_t x = mont::sub(mont::sub(mont::square(slope), p.x), q.x);
    uint131_t y = mont::sub(mont::mul(slope, mont::sub(p.x, x)), p.y);

    _x[i] = x;
    _y[i] = y;
    _walk_len[i]++;
  }
}

void CPUPointFinderPrime::init()
{
  init("");
}

void CPUPointFinderPrime::init(const std::string& file)
{
  load(file);
  _chain.resize(_num_points);
}

double CPUPointFinderPrime::step()
{
  auto start = std::chrono::steady_clock::now();
  std::vector<DistinguishedPoint> results;

  for(int iteration = 0; iteration < _iters_per_step; iteration++) {
    do_step_prime(results);
    _counter++;
  }

  if(!results.empty() && _callback) {
    _callback(results);
  }

  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
