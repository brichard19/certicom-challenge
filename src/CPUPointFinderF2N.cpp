#include "CPUPointFinderF2N.h"

#include <cassert>
#include <chrono>
#include <vector>

#include "ecc_internal.h"
#include "gf2.h"

CPUPointFinderF2N::CPUPointFinderF2N(int dpbits, size_t num_points, bool benchmark)
    : CPUPointFinder(dpbits, num_points), _benchmark(benchmark)
{
}

void CPUPointFinderF2N::init() { init(""); }

void CPUPointFinderF2N::init(const std::string& file)
{
  load(file);
  _chain.resize(_num_points);
}

void CPUPointFinderF2N::do_step_binary(std::vector<DistinguishedPoint>& results)
{
  const uint131_t one = make_uint131(1);

  for(size_t i = 0; i < _num_points; i++) {
    if(_benchmark == false && (_x[i].w.v0 & _dpmask) == 0) {
      ecc::ecpoint_t point(_x[i], _y[i]);
      assert(ecc::exists(point));
      results.emplace_back(_priv[i], point, _dpbits, _walk_len[i]);
      generate_walk(i);
    }

    size_t rw_index = _x[i].w.v0 & 0x1f;
    uint131_t denominator = gf2::add(_rx[rw_index], _x[i]);
    _chain[i] = gf2::mul(i == 0 ? one : _chain[i - 1], denominator, _params.field);
  }

  uint131_t inverse = gf2::inv(_chain.empty() ? one : _chain.back(), _params.field);

  for(size_t i = _num_points; i-- > 0;) {
    size_t rw_index = _x[i].w.v0 & 0x1f;
    uint131_t px = _x[i];
    uint131_t py = _y[i];
    uint131_t qx = _rx[rw_index];
    uint131_t qy = _ry[rw_index];

    uint131_t denominator = gf2::add(qx, px);
    uint131_t denominator_inverse =
        i == 0 ? inverse : gf2::mul(inverse, _chain[i - 1], _params.field);
    inverse = gf2::mul(inverse, denominator, _params.field);

    uint131_t lambda = gf2::mul(gf2::add(qy, py), denominator_inverse, _params.field);
    uint131_t x = gf2::add(gf2::square(lambda, _params.field), lambda);
    x = gf2::add(gf2::add(gf2::add(x, px), qx), _params.a);
    uint131_t y = gf2::mul(lambda, gf2::add(px, x), _params.field);
    y = gf2::add(gf2::add(y, x), py);

    _x[i] = x;
    _y[i] = y;
    _walk_len[i]++;
  }
}

double CPUPointFinderF2N::step()
{
  auto start = std::chrono::steady_clock::now();
  std::vector<DistinguishedPoint> results;

  for(int iteration = 0; iteration < _iters_per_step; iteration++) {
    do_step_binary(results);
    _counter++;
  }

  if(!results.empty() && _callback) {
    _callback(results);
  }

  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
