#include "CPUPointFinder.h"
#include "CPUPointFinderF2N.h"

#include <cassert>
#include <chrono>
#include <fstream>
#include <stdexcept>

#include "ecc_internal.h"
#include "gf2.h"
#include "montgomery.h"
#include "util.h"

namespace {

struct VecUint131 {
  uint8_t data[17];
};

void store_uint131(void* ptr, size_t index, size_t count, const uint131_t& value)
{
  uint128_t* ptr128 = static_cast<uint128_t*>(ptr);
  uint8_t* ptr8 = static_cast<uint8_t*>(ptr);
  ptr128[index] = (uint128_t(value.w.v1) << 64) | value.w.v0;
  ptr8[count * sizeof(uint128_t) + index] = uint8_t(value.w.v2);
}

uint131_t load_uint131(const void* ptr, size_t index, size_t count)
{
  const uint128_t* ptr128 = static_cast<const uint128_t*>(ptr);
  const uint8_t* ptr8 = static_cast<const uint8_t*>(ptr);
  uint128_t low = ptr128[index];
  uint131_t value = {};
  value.w.v0 = uint64_t(low);
  value.w.v1 = uint64_t(low >> 64);
  value.w.v2 = ptr8[count * sizeof(uint128_t) + index];
  return value;
}

} // namespace

// Common helpers implemented using the protected members of CPUPointFinder

void checked_read(std::ifstream& file, void* ptr, size_t size)
{
  if(!file.read(static_cast<char*>(ptr), size)) {
    throw std::runtime_error("Error reading point-finder progress file");
  }
}

CPUPointFinder::CPUPointFinder(int dpbits, size_t num_points)
    : _dpbits(dpbits), _num_points(num_points)
{
  if(dpbits < 1 || dpbits > 63) {
    throw std::invalid_argument("dpbits must be between 1 and 63");
  }
  if(num_points == 0) {
    throw std::invalid_argument("num_points must be greater than zero");
  }
  _dpmask = (uint64_t(1) << dpbits) - 1;
}

void CPUPointFinder::generate_walk(size_t index)
{
  ecc::ecpoint_t point;
  do {
    _priv[index] = ecc::genkey();
    point = ecc::mul(_priv[index], ecc::g());
  } while(ecc::is_infinity(point));
  _x[index] = point.x;
  _y[index] = point.y;
  _walk_len[index] = 0;
}

void CPUPointFinder::load(const std::string& file)
{
  std::vector<RWPoint> rw = get_rw_points();
  _rx.resize(rw.size());
  _ry.resize(rw.size());
  for(size_t i = 0; i < rw.size(); i++) {
    _rx[i] = rw[i].p.x;
    _ry[i] = rw[i].p.y;
  }

  if(!file.empty() && util::file_exists(file)) {
    // load progress
    std::ifstream f(file, std::ios::binary);
    if(!f)
      throw std::runtime_error("Unable to open progress file for reading");

    checked_read(f, &_num_points, sizeof(_num_points));
    checked_read(f, &_counter, sizeof(_counter));

    _x.resize(_num_points);
    _y.resize(_num_points);
    _priv.resize(_num_points);
    _walk_len.resize(_num_points);

    checked_read(f, _priv.data(), sizeof(uint131_t) * _num_points);

    std::vector<uint8_t> packed(sizeof(VecUint131) * _num_points);
    checked_read(f, packed.data(), packed.size());
    for(size_t i = 0; i < _num_points; i++)
      _x[i] = load_uint131(packed.data(), i, _num_points);

    checked_read(f, packed.data(), packed.size());
    for(size_t i = 0; i < _num_points; i++)
      _y[i] = load_uint131(packed.data(), i, _num_points);

    std::vector<uint64_t> walk_start(_num_points);
    checked_read(f, walk_start.data(), sizeof(uint64_t) * walk_start.size());
    for(size_t i = 0; i < _num_points; i++) {
      if(walk_start[i] > _counter) {
        throw std::runtime_error("Invalid walk start in progress file");
      }
      _walk_len[i] = _counter - walk_start[i];
      if(!ecc::exists(ecc::ecpoint_t(_x[i], _y[i]))) {
        throw std::runtime_error("Invalid point in progress file");
      }
    }

    return;
  }

  _x.resize(_num_points);
  _y.resize(_num_points);
  _priv.resize(_num_points);
  _walk_len.resize(_num_points);

  for(size_t i = 0; i < _num_points; i++) {
    _priv[i] = ecc::genkey();
  }

  std::vector<ecc::ecpoint_t> points = ecc::mul(_priv, ecc::g());
  for(size_t i = 0; i < _num_points; i++) {
    if(ecc::is_infinity(points[i])) {
      generate_walk(i);
    } else {
      _x[i] = points[i].x;
      _y[i] = points[i].y;
      _walk_len[i] = 0;
    }
  }
}

void CPUPointFinder::save_progress_impl(const std::string& file_name)
{
  std::ofstream file(file_name, std::ios::binary);
  if(!file) {
    throw std::runtime_error("Unable to open progress file for writing");
  }

  file.write(reinterpret_cast<const char*>(&_num_points), sizeof(_num_points));
  file.write(reinterpret_cast<const char*>(&_counter), sizeof(_counter));
  file.write(reinterpret_cast<const char*>(_priv.data()), sizeof(uint131_t) * _num_points);

  std::vector<uint8_t> packed(sizeof(VecUint131) * _num_points);
  for(size_t i = 0; i < _num_points; i++)
    store_uint131(packed.data(), i, _num_points, _x[i]);
  file.write(reinterpret_cast<const char*>(packed.data()), packed.size());

  for(size_t i = 0; i < _num_points; i++)
    store_uint131(packed.data(), i, _num_points, _y[i]);
  file.write(reinterpret_cast<const char*>(packed.data()), packed.size());

  std::vector<uint64_t> walk_start(_num_points);
  for(size_t i = 0; i < _num_points; i++)
    walk_start[i] = _counter - _walk_len[i];
  file.write(reinterpret_cast<const char*>(walk_start.data()),
             sizeof(uint64_t) * walk_start.size());

  if(!file) {
    throw std::runtime_error("Error writing point-finder progress file");
  }
}

void CPUPointFinder::set_callback(
    std::function<void(const std::vector<DistinguishedPoint>&)> callback)
{
  _callback = callback;
}

void CPUPointFinder::save_progress(const std::string& file) { save_progress_impl(file); }

size_t CPUPointFinder::work_per_step() { return _num_points * _iters_per_step; }

int CPUPointFinder::iters_per_step() { return _iters_per_step; }

int CPUPointFinder::parallel_walks() { return int(_num_points); }

// Base class is abstract; no default init/step implementations here.

// Define destructor to emit vtable in this translation unit.
CPUPointFinder::~CPUPointFinder() = default;

// Factory implementation
std::unique_ptr<CPUPointFinder> make_cpu_point_finder(int dpbits, size_t num_points)
{
  return std::make_unique<CPUPointFinderF2N>(dpbits, num_points);
}
