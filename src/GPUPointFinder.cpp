#include "fmt/format.h"
#include <assert.h>
#include <chrono>
#include <fstream>
#include <hip/hip_runtime.h>
#include <map>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <stdint.h>

#include "GPUPointFinder.h"
#include "ec_rho.h"
#include "ecc.cuh"
#include "hip_helper.h"
#include "log.h"
#include "montgomery.h"
#include "util.h"

namespace {

#define IFSTREAM_CALL(condition)                                                                   \
  {                                                                                                \
    if(!condition) {                                                                               \
      std::stringstream ss;                                                                        \
      ss << "FILE error " << __LINE__ << std::endl;                                                \
      throw std::runtime_error(ss.str());                                                          \
    }                                                                                              \
  }

struct KernelInfo {
  void *do_step;
  void *batch_mul;
  void *sanity_check;
  void *refill_staging;
};

std::map<std::string, KernelInfo> _kernel_map = {
    {"ecp79", {(void *)do_step_p79, (void *)batch_multiply_p79, (void *)sanity_check_p79}},
    {"ecp89", {(void *)do_step_p89, (void *)batch_multiply_p89, (void *)sanity_check_p89}},
    {"ecp131", {(void *)do_step_p131, (void *)batch_multiply_p131, (void *)sanity_check_p131}},
};

uint131_t load_uint131(const void *p, int idx, int n)
{
  const uint128_t *p128 = (const uint128_t *)p;
  const uint8_t *p8 = (const uint8_t *)p;

  uint128_t u128 = p128[idx];

  uint131_t x;
  x.w.v0 = (uint64_t)u128;
  x.w.v1 = (uint64_t)(u128 >> 64);
  x.w.v2 = (uint32_t)p8[n * sizeof(uint128_t) + idx];

  return x;
}

void store_uint131(void *p, int idx, int n, uint131_t x)
{
  uint128_t *p128 = (uint128_t *)p;
  uint8_t *p8 = (uint8_t *)p;

  uint128_t u128 = ((uint128_t)x.w.v1 << 64) | x.w.v0;
  p128[idx] = u128;

  p8[n * sizeof(uint128_t) + idx] = (uint8_t)x.w.v2;
}

} // namespace

// uint131_t is 20 bytes while we only need 17 bytes to store it, so
// to allocate vectorized uint131_t arrays, use sizeof(vec_uint131_t);
typedef struct {
  uint8_t data[17];
} vec_uint131_t;

GPUPointFinder::GPUPointFinder(int device, int dpbits, bool benchmark)
{
  assert(dpbits >= 16 && dpbits <= 63);

  _benchmark = benchmark;

  std::string curve_name = ecc::curve_name();

  _do_step_ptr = _kernel_map[curve_name].do_step;
  _batch_multiply_ptr = _kernel_map[curve_name].batch_mul;
  _sanity_check_ptr = _kernel_map[curve_name].sanity_check;

  _device = device;
  _dpbits = dpbits;

  int compute_units = get_cu_count(device);
  int simd_width = get_warp_size(device);

  // Use these for now
  _blocks = compute_units * 32;
  _threads_per_block = simd_width;

  HIP_CALL(hipSetDevice(device));
  HIP_CALL(hipSetDeviceFlags(hipDeviceScheduleYield));

#ifdef DEBUG
  LOG("Debug mode ON");
#endif

  LOG("Initializing GPU {}: {}", _device, get_gpu_name(_device));
  LOG("Compute units: {}", compute_units);
  LOG("SIMD width:    {}", simd_width);

  int num_points = _blocks * _threads_per_block * POINTS_PER_THREAD;

  uint32_t total_threads = _blocks * _threads_per_block;

  // Round up to multiple of number of threads to ensure CUs get roughtly
  // the same amount of work
  _num_points = (num_points + total_threads - 1) / total_threads * total_threads;

  // Create DP mask
  _dpmask = ((uint64_t)1 << _dpbits) - 1;

  // Using the number of DP bits and the number of parallel walks, we can
  // approximate how many results we will get per iteration
  double prob = 1.0 / pow(2.0, (double)_dpbits);

  const int iters = 64 * 1024;
  // Can run ~65k iteration before buffer is full
  _result_buf_size = iters * (int)((double)_num_points * prob + 1.0);

  _report_count = (int)((double)_result_buf_size * 0.85);

  // Calculate an appropriate size for the staging buffer
  int staging_buf_size = (int)(((double)_num_points * prob + 1.0) * iters);

  _staging = stack_create<StagingPoint>(staging_buf_size);
}

GPUPointFinder::~GPUPointFinder()
{
  stack_destroy(_staging);
  free_buffers();
}

void GPUPointFinder::report_points()
{
  int count = *_result_count;

  if(count > 0) {

    std::vector<DistinguishedPoint> dps;
    for(int i = 0; i < count; i++) {
      uint131_t x = _result_buf[i].x;
      uint131_t y = _result_buf[i].y;
      uint131_t a = _result_buf[i].a;
      uint64_t length = _result_buf[i].length;

      ecc::ecpoint_t p(x, y);

      // Do a quick verification
      assert(ecc::exists(p));

      assert((p.x.w.v0 & _dpmask) == 0);

      dps.push_back(DistinguishedPoint(a, p, _dpbits, length));
    }

    if(_callback) {
      _callback(dps);
    }

    *_result_count = 0;
  }
}

void GPUPointFinder::refill_staging()
{
  LOG("Refilling staging buffer");

  int n = _staging.max_size - *_staging.size;

  // Create lookup table for G, 2G, 4G ... nG
  uint131_t *dev_gx = nullptr;
  uint131_t *dev_gy = nullptr;
  uint131_t *mbuf = nullptr;
  uint131_t *priv_keys = nullptr;
  uint131_t *x_ptr = nullptr;
  uint131_t *y_ptr = nullptr;

  HIP_CALL(hipMallocManaged(&dev_gx, sizeof(uint131_t) * (ecc::curve_strength() + 1)));
  HIP_CALL(hipMallocManaged(&dev_gy, sizeof(uint131_t) * (ecc::curve_strength() + 1)));
  HIP_CALL(hipMalloc(&mbuf, sizeof(uint131_t) * n));

  HIP_CALL(hipMallocManaged(&priv_keys, sizeof(uint131_t) * n));
  HIP_CALL(hipMallocManaged(&x_ptr, sizeof(uint131_t) * n));
  HIP_CALL(hipMallocManaged(&y_ptr, sizeof(uint131_t) * n));

  // Generate G, 2G, 4G, 8G ... nG for batch multiplication
  ecc::ecpoint_t g = ecc::g();

  for(int i = 0; i < ecc::curve_strength() + 1; i++) {
    dev_gx[i] = g.x;
    dev_gy[i] = g.y;

    g = ecc::dbl(g);
  }

  // Generate random keys
  for(int i = 0; i < n; i++) {
    priv_keys[i] = ecc::genkey();
  }

  // Initialize keys
  HIP_CALL(hipLaunchKernel((void *)clear_public_keys, dim3(_blocks), dim3(_threads_per_block), 0,
                           x_ptr, y_ptr, n));

  // Initialize keys
  int bits = ecc::curve_strength();
  for(int b = 0; b < bits; b++) {
    HIP_CALL(hipLaunchKernel(_batch_multiply_ptr, dim3(_blocks), dim3(_threads_per_block), 0, x_ptr,
                             y_ptr, priv_keys, mbuf, dev_gx, dev_gy, b, n));
  }
  HIP_CALL(hipDeviceSynchronize());

  if(_verify_points) {
    // Verify
    LOG("Verifying");
    std::vector<uint131_t> keys;
    for(int i = 0; i < n; i++) {
      keys.push_back(priv_keys[i]);
    }

    auto expected = ecc::mul(keys, ecc::g());

    for(int i = 0; i < n; i++) {
      ecc::ecpoint_t actual = ecc::ecpoint_t(load_uint131(x_ptr, i, n), load_uint131(y_ptr, i, n));

      assert(ecc::is_equal(expected[i], actual));
    }
  }

  std::vector<StagingPoint> tmp;
  // Copy into staging buffer
  for(int i = 0; i < n; i++) {
    StagingPoint sp;
    sp.a = priv_keys[i];
    sp.x = load_uint131(x_ptr, i, n);
    sp.y = load_uint131(y_ptr, i, n);
    tmp.push_back(sp);
  }

  stack_fill(_staging, tmp.data(), tmp.size());

  HIP_CALL(hipFree(dev_gx));
  HIP_CALL(hipFree(dev_gy));
  HIP_CALL(hipFree(mbuf));
  HIP_CALL(hipFree(x_ptr));
  HIP_CALL(hipFree(y_ptr));
  HIP_CALL(hipFree(priv_keys));
}

void GPUPointFinder::init() { init(""); }

void GPUPointFinder::init(const std::string &filename)
{
  allocate_buffers(_num_points);

  // Don't need staging buffer for benchmark
  if(_benchmark == false) {
    refill_staging();
  }

  // Initialize random walk points
  std::vector<RWPoint> rw = get_rw_points();

  std::vector<uint8_t> rx(sizeof(vec_uint131_t) * rw.size());
  std::vector<uint8_t> ry(sizeof(vec_uint131_t) * rw.size());

  for(int i = 0; i < rw.size(); i++) {
    store_uint131(rx.data(), i, rw.size(), rw[i].p.x);
    store_uint131(ry.data(), i, rw.size(), rw[i].p.y);
  }

  HIP_CALL(hipMemcpy(_dev_rx, rx.data(), rx.size() * sizeof(rx[0]), hipMemcpyHostToDevice));
  HIP_CALL(hipMemcpy(_dev_ry, ry.data(), ry.size() * sizeof(ry[0]), hipMemcpyHostToDevice));

  if(!filename.empty() && util::file_exists(filename)) {
    LOG("Loading {}", filename);
    load(filename);

    // Check sanity
    *_sanity_flag = 0;

    HIP_CALL(hipLaunchKernel((void *)_sanity_check_ptr, dim3(_blocks), dim3(_threads_per_block), 0,
                             _dev_x, _dev_y, _num_points, _sanity_flag));
    HIP_CALL(hipDeviceSynchronize());

    if(*_sanity_flag != 0) {
      LOG("Detected {} errors", *_sanity_flag);
      throw std::runtime_error("Error verifying points");
    }

  } else {

    if(!filename.empty()) {
      LOG("Data will be saved to: {}", filename);
    }

    // Generate new private keys
    for(int i = 0; i < _num_points; i++) {
      _priv_key_a[i] = ecc::genkey();
    }

    // Clear keys on device
    HIP_CALL(hipLaunchKernel((void *)clear_public_keys, dim3(_blocks), dim3(_threads_per_block), 0,
                             _dev_x, _dev_y, _num_points));
    HIP_CALL(hipDeviceSynchronize());

    HIP_CALL(hipLaunchKernel((void *)reset_counters, dim3(_blocks), dim3(_threads_per_block), 0,
                             _walk_start, _counter, _num_points));
    HIP_CALL(hipDeviceSynchronize());

    uint131_t *dev_gx = nullptr;
    uint131_t *dev_gy = nullptr;

    HIP_CALL(hipMallocManaged(&dev_gx, sizeof(uint131_t) * (ecc::curve_strength() + 1)));
    HIP_CALL(hipMallocManaged(&dev_gy, sizeof(uint131_t) * (ecc::curve_strength() + 1)));

    // Generate G, 2G, 4G, 8G ... nG for batch multiplication
    ecc::ecpoint_t g = ecc::g();

    for(int i = 0; i < ecc::curve_strength() + 1; i++) {
      dev_gx[i] = g.x;
      dev_gy[i] = g.y;

      g = ecc::dbl(g);
    }

    // Initialize keys
    int bits = ecc::curve_strength();
    for(int i = 0; i < bits; i++) {
      HIP_CALL(hipLaunchKernel(_batch_multiply_ptr, dim3(_blocks), dim3(_threads_per_block), 0,
                               _dev_x, _dev_y, _priv_key_a, _mbuf, dev_gx, dev_gy, i, _num_points));
      HIP_CALL(hipDeviceSynchronize());
    }

    HIP_IGNORE(hipFree(dev_gx));
    HIP_IGNORE(hipFree(dev_gy));

    // Verify
    if(_verify_points) {
      std::vector<uint8_t> x(_num_points * sizeof(vec_uint131_t));
      std::vector<uint8_t> y(_num_points * sizeof(vec_uint131_t));

      HIP_CALL(hipMemcpy(x.data(), _dev_x, x.size(), hipMemcpyDeviceToHost));
      HIP_CALL(hipMemcpy(y.data(), _dev_y, y.size(), hipMemcpyDeviceToHost));

      for(int i = 0; i < 1000; i++) {
        uint131_t k1 = _priv_key_a[i];

        ecc::ecpoint_t sum = ecc::mul(k1, ecc::g());

        uint131_t xval = load_uint131(x.data(), i, _num_points);
        uint131_t yval = load_uint131(y.data(), i, _num_points);

        if(!ecc::exists(ecc::ecpoint_t(xval, yval))) {
          LOG("{}", to_str(sum.x));
          LOG("{}", to_str(xval));
          LOG("Validation failed: Point does not exist");
          exit(1);
        }

        if(!ecc::is_equal(sum, ecc::ecpoint_t(xval, yval))) {
          LOG("Expected:");
          LOG("{}", to_str(sum.x));
          LOG("{}", to_str(sum.y));
          LOG("{}", to_str(mont::from(sum.x)));
          LOG("{}", to_str(mont::from(sum.y)));

          LOG("Got:");
          LOG("{}", to_str(xval));
          LOG("{}", to_str(yval));
          LOG("{}", to_str(mont::from(xval)));
          LOG("{}", to_str(mont::from(yval)));

          LOG("Validation failed: Point exists but is incorrect for key");
          exit(1);
        }
      }
    }
  }
}

void GPUPointFinder::free_buffers()
{
  HIP_IGNORE(hipFree(_mbuf));
  HIP_IGNORE(hipFree(_dev_x));
  HIP_IGNORE(hipFree(_dev_y));
  HIP_IGNORE(hipFree(_dev_rx));
  HIP_IGNORE(hipFree(_dev_ry));

  HIP_IGNORE(hipFree(_priv_key_a));
  HIP_IGNORE(hipFree(_result_buf));
  HIP_IGNORE(hipFree(_result_count));
  HIP_IGNORE(hipFree(_sanity_flag));

  HIP_IGNORE(hipFree(_walk_start));
}

void GPUPointFinder::allocate_buffers(int n)
{
  free_buffers();

  HIP_CALL(hipMalloc(&_mbuf, n * sizeof(vec_uint131_t)));
  HIP_CALL(hipMalloc(&_dev_x, n * sizeof(vec_uint131_t)));
  HIP_CALL(hipMalloc(&_dev_y, n * sizeof(vec_uint131_t)));
  HIP_CALL(hipMalloc(&_dev_rx, _NUM_R_POINTS * sizeof(vec_uint131_t)));
  HIP_CALL(hipMalloc(&_dev_ry, _NUM_R_POINTS * sizeof(vec_uint131_t)));

  HIP_CALL(hipMallocManaged(&_priv_key_a, _num_points * sizeof(uint131_t)));

  HIP_CALL(hipMallocManaged(&_result_buf, _result_buf_size * sizeof(DPResult)));

  HIP_CALL(hipMallocManaged(&_result_count, sizeof(uint32_t)));
  *_result_count = 0;

  HIP_CALL(hipMallocManaged(&_sanity_flag, sizeof(uint32_t)));
  *_sanity_flag = 0;

  HIP_CALL(hipMallocManaged(&_walk_start, _num_points * sizeof(uint64_t)));
}

double GPUPointFinder::step()
{

  // We don't need result buf for benchmark
  DPResult *result_buf = _benchmark ? nullptr : _result_buf;

  hipEvent_t start;
  hipEvent_t stop;
  float elapsed = 0.0f;

  HIP_CALL(hipEventCreate(&start));
  HIP_CALL(hipEventCreate(&stop));

  HIP_CALL(hipEventRecord(start));
  for(int i = 0; i < _iters_per_step; i++) {
    HIP_CALL(hipLaunchKernel(_do_step_ptr, dim3(_blocks), dim3(_threads_per_block), 0, _dev_x,
                             _dev_y, _dev_rx, _dev_ry, _mbuf, _num_points, result_buf,
                             _result_count, _staging, _priv_key_a, _counter, _walk_start, _dpmask));
    _counter++;
  }
  HIP_CALL(hipEventRecord(stop));

  HIP_CALL(hipDeviceSynchronize());

  HIP_CALL(hipEventElapsedTime(&elapsed, start, stop));

  if(_first_run || _verify_points) {
    HIP_CALL(hipLaunchKernel((void *)_sanity_check_ptr, dim3(_blocks), dim3(_threads_per_block), 0,
                             _dev_x, _dev_y, _num_points, _sanity_flag));
    HIP_CALL(hipDeviceSynchronize());

    if(*_sanity_flag != 0) {
      LOG("Detected {} errors\n", *_sanity_flag);

      std::vector<uint8_t> x(_num_points * sizeof(vec_uint131_t));
      std::vector<uint8_t> y(_num_points * sizeof(vec_uint131_t));

      HIP_CALL(hipMemcpy(x.data(), _dev_x, x.size(), hipMemcpyDeviceToHost));
      HIP_CALL(hipMemcpy(y.data(), _dev_y, y.size(), hipMemcpyDeviceToHost));

      int count = 0;
      for(int i = 0; i < _num_points; i++) {
        uint131_t xval = load_uint131(x.data(), i, _num_points);
        uint131_t yval = load_uint131(y.data(), i, _num_points);
        ecc::ecpoint_t p(xval, yval);

        if(!ecc::exists(p)) {
          count++;
          LOG("{} {}", to_str(p.x), to_str(p.y));
        }
      }

      LOG("Host detected {} errors", count);

      throw std::runtime_error("Error verifying points");
    }

    // Only verify on the first iteration
    _first_run = false;
  }

  // Don't need to check these during benchmark
  if(_benchmark == false) {

    // Refill at 10% capacity
    if(*_staging.size <= (int)(_staging.max_size * 0.1)) {
      refill_staging();
    }

    if(*_result_count >= _report_count) {
      report_points();
    }
  }

  return (double)elapsed / 1000.0f;
}

void GPUPointFinder::set_callback(
    std::function<void(const std::vector<DistinguishedPoint> &)> callback)
{
  _callback = callback;
}

size_t GPUPointFinder::work_per_step() { return _num_points * _iters_per_step; }

int GPUPointFinder::iters_per_step() { return _iters_per_step; }

int GPUPointFinder::parallel_walks() { return _num_points; }

void GPUPointFinder::save_progress(const std::string &file_name)
{

  report_points();

  LOG("Saving progress to {}", file_name);

  std::ofstream file;

  file.open(file_name, std::ios::binary);

  if(!file.good()) {
    LOG("File not good");
    return;
  }

  size_t count = _num_points;

  file.write((char *)&count, sizeof(count));
  if(!file.good()) {
    LOG("File not good");
    return;
  }

  // Write the counter
  file.write((char *)&_counter, sizeof(_counter));
  if(!file.good()) {
    LOG("File not good");
    return;
  }

  // Private keys a
  file.write((char *)_priv_key_a, sizeof(uint131_t) * count);
  if(!file.good()) {
    LOG("File not good");
    return;
  }

  std::vector<char> tmp(sizeof(vec_uint131_t) * count);
  // Points X
  HIP_CALL(hipMemcpy(tmp.data(), _dev_x, sizeof(vec_uint131_t) * count, hipMemcpyDeviceToHost));
  file.write(tmp.data(), sizeof(vec_uint131_t) * count);

  // Points Y
  HIP_CALL(hipMemcpy(tmp.data(), _dev_y, sizeof(vec_uint131_t) * count, hipMemcpyDeviceToHost));
  file.write(tmp.data(), sizeof(vec_uint131_t) * count);

  // Write starting counters
  memcpy(tmp.data(), _walk_start, sizeof(uint64_t) * count);
  file.write(tmp.data(), sizeof(uint64_t) * count);

  file.close();
}

void GPUPointFinder::load(const std::string &file_name)
{
  std::ifstream file;

  file.open(file_name, std::ios::in | std::ios::binary);

  if(!file.is_open()) {
    throw std::runtime_error("Unable to open file for reading");
  }

  // Load the count
  size_t count = 0;

  IFSTREAM_CALL(file.read((char *)&count, sizeof(count)));

  assert(count == _num_points);

  // Load counter
  IFSTREAM_CALL(file.read((char *)&_counter, sizeof(_counter)));

  IFSTREAM_CALL(file.read((char *)_priv_key_a, sizeof(uint131_t) * count));

  std::vector<char> tmp(sizeof(vec_uint131_t) * count);

  IFSTREAM_CALL(file.read(tmp.data(), sizeof(vec_uint131_t) * count));
  HIP_CALL(hipMemcpy(_dev_x, tmp.data(), sizeof(vec_uint131_t) * count, hipMemcpyHostToDevice));

  IFSTREAM_CALL(file.read(tmp.data(), sizeof(vec_uint131_t) * count));
  HIP_CALL(hipMemcpy(_dev_y, tmp.data(), sizeof(vec_uint131_t) * count, hipMemcpyHostToDevice));

  IFSTREAM_CALL(file.read(tmp.data(), sizeof(uint64_t) * count));
  memcpy(_walk_start, tmp.data(), sizeof(uint64_t) * count);
}