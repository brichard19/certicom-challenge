#include <assert.h>
#include <chrono>
#include <fstream>
#include <hip/hip_runtime.h>
#include <sstream>
#include <stdexcept>
#include <stdint.h>

#include "GPUPointFinder.h"
#include "ec_rho.h"
#include "ecc.cuh"
#include "hip_helper.h"
#include "log.h"
#include "mont.h"
#include "util.h"

namespace {

#define IFSTREAM_CALL(condition)\
{\
  if(!condition) {\
    std::stringstream ss;\
    ss << "FILE error " << __LINE__ << std::endl;\
    throw std::runtime_error(ss.str());\
  }\
}\

template<class... T> hipError_t hipLaunchKernel(void* kernel, dim3 gridDim, dim3 blockDim, size_t sharedMem, T... args)
{
  std::vector<void*> ptr = {&args...};

  return hipLaunchKernel(kernel, gridDim, blockDim, ptr.data(), sharedMem, 0);
}

std::string get_gpu_name(int device_id)
{
  hipDeviceProp_t props;
 
  HIP_CALL(hipGetDeviceProperties(&props, device_id));

  return std::string(props.name);
}

int get_cu_count(int device_id)
{
  hipDeviceProp_t props;
 
  HIP_CALL(hipGetDeviceProperties(&props, device_id));

  return props.multiProcessorCount;
}

int get_warp_size(int device_id)
{
  hipDeviceProp_t props;
 
  HIP_CALL(hipGetDeviceProperties(&props, device_id));

  return props.warpSize;
}

}

GPUPointFinder::GPUPointFinder(int device, uint32_t num_points, int dpbits, bool benchmark)
{
  assert(dpbits >= 16 && dpbits <= 63);

  _benchmark = benchmark;

  std::string curve_name = ecc::curve_name();

  if(curve_name == "ecp131") {
    _do_step_ptr = (void*)do_step_p131;
    _batch_multiply_ptr = (void*)batch_multiply_p131;
    _sanity_check_ptr = (void*)sanity_check_p131;
  } else if(curve_name == "ecp79") {
    _do_step_ptr = (void*)do_step_p79;
    _batch_multiply_ptr = (void*)batch_multiply_p79;
    _sanity_check_ptr = (void*)sanity_check_p79;
  } else {
    LOG("Invalid curve!");
    throw std::runtime_error("Invalid curve");
  }

  _device = device;
  _dpbits = dpbits;

  int compute_units = get_cu_count(device);
  int simd_width = get_warp_size(device);

  // Use these for now
  _blocks = compute_units * 8;
  _threads = simd_width;
  
  HIP_CALL(hipSetDevice(device));
  LOG("Initializing GPU {}: {}", _device, get_gpu_name(_device));
  LOG("Compute units: {}", compute_units);
  LOG("SIMD width:    {}", simd_width);

  // Give 256 points to each thread unless specified
  if(num_points == 0) {
    num_points = _blocks * _threads * 256;
  }

  uint32_t total_threads = _blocks * _threads;

  // Round up to multiple of number of threads to ensure CUs get roughtly
  // the same amount of work
  _num_points = (num_points + total_threads - 1) / total_threads * total_threads;

  // Create DP mask
  _dpmask = ((uint64_t)1 << _dpbits) - 1;

  // Using the number of DP bits and the number of parallel walks, we can
  // approximate how many results we will get per iteration
  double prob = 1.0 / pow(2, _dpbits);

  // Can run ~65k iteration before buffer is full
  _result_buf_size = 65536 * (int)((double)_num_points * prob + 1.0);

  _report_count = (int)((double)_result_buf_size * 0.85);

  // Calculate an appropriate size for the staging buffer
  _staging_buf_size = (int)(((double)_num_points * prob + 1.0) * 65536);

  // Refill after we reach 15%
  _staging_min = (int)((double)_staging_buf_size * 0.15);
}


GPUPointFinder::~GPUPointFinder()
{
  free_buffers();
}

void GPUPointFinder::load(const std::string& file_name)
{
  std::ifstream file;

  file.open(file_name, std::ios::in | std::ios::binary);

  if(!file.is_open()) {
    throw std::runtime_error("Unable to open file for reading");
  }

  // Load the count
  size_t count = 0;

  IFSTREAM_CALL(file.read((char*)&count, sizeof(count)));
  
  assert(count == _num_points);

  // Load counter
  IFSTREAM_CALL(file.read((char*)&_counter, sizeof(_counter)));

  IFSTREAM_CALL(file.read((char*)_priv_key_a, sizeof(uint131_t) * count));

  std::vector<char> tmp(sizeof(uint131_t) * count);

  IFSTREAM_CALL(file.read(tmp.data(), sizeof(uint131_t) * count));
  HIP_CALL(hipMemcpy(_dev_x, tmp.data(), sizeof(uint131_t) * count, hipMemcpyHostToDevice));
  
  IFSTREAM_CALL(file.read(tmp.data(), sizeof(uint131_t) * count));
  HIP_CALL(hipMemcpy(_dev_y, tmp.data(), sizeof(uint131_t) * count, hipMemcpyHostToDevice));
  
  IFSTREAM_CALL(file.read(tmp.data(), sizeof(uint64_t) * count));
  memcpy(_walk_start, tmp.data(), sizeof(uint64_t) * count);

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
 
      dps.push_back(DistinguishedPoint(a, p, length));
    }

    if(_callback) {
      _callback(dps);
    }

    *_result_count = 0;
  }
}

void GPUPointFinder::refill_staging()
{
  int count = *_staging_count;

  LOG("Refilling staging buffer");

  for(int i = count; i < _staging_buf_size; i++) {
    StagingPoint p;

    ecc::ecpoint_t point;

    // Generate points, avoid distinguished points
    do {
      p.a = ecc::genkey();

      point = ecc::mul(p.a, ecc::g());
    } while((point.x.w.v0 & _dpmask) == 0);

    assert(ecc::exists(point));
    
    p.x = point.x;
    p.y = point.y;

    _staging_buf[i] = p;
  }

  *_staging_count= _staging_buf_size;
}

void GPUPointFinder::init()
{
  init("");
}

void GPUPointFinder::init(const std::string& filename)
{  
  allocate_buffers(_num_points);

  // Don't need staging for benchmark, avoid costly refill
  if(_benchmark == false) {
    refill_staging();
  }

  // Initialize random walk points
  std::vector<RWPoint> rw = get_rw_points();

  std::vector<uint131_t> rx;
  std::vector<uint131_t> ry;
  for(auto p : rw) {
    rx.push_back(p.p.x);
    ry.push_back(p.p.y);
  }

  HIP_CALL(hipMemcpy(_dev_rx, rx.data(), rx.size() * sizeof(rx[0]), hipMemcpyHostToDevice));
  HIP_CALL(hipMemcpy(_dev_ry, ry.data(), ry.size() * sizeof(ry[0]), hipMemcpyHostToDevice));


  if(!filename.empty() && util::file_exists(filename)) {
    LOG("Loading {}", filename);
    load(filename);

    // Check sanity
    *_sanity_flag = 0;

    HIP_CALL(hipLaunchKernel((void*)_sanity_check_ptr, dim3(_blocks), dim3(_threads), 0, _dev_x, _dev_y, _num_points, _sanity_flag));
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
    HIP_CALL(hipLaunchKernel((void*)clear_public_keys, dim3(_blocks), dim3(_threads), 0, _dev_x, _dev_y, _num_points));
    HIP_CALL(hipDeviceSynchronize());

    HIP_CALL(hipLaunchKernel((void*)reset_counters, dim3(_blocks), dim3(_threads), 0, _walk_start, _counter, _num_points));
    HIP_CALL(hipDeviceSynchronize());

    uint131_t* dev_gx = nullptr;
    uint131_t* dev_gy = nullptr;

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
      HIP_CALL(hipLaunchKernel(_batch_multiply_ptr, dim3(_blocks), dim3(_threads), 0, _dev_x, _dev_y, _priv_key_a, _mbuf, dev_gx, dev_gy, i, _num_points));
      HIP_CALL(hipDeviceSynchronize());
    }

    HIP_IGNORE(hipFree(dev_gx));
    HIP_IGNORE(hipFree(dev_gy));

    // Verify
    if(_verify_points) {
      std::vector<uint131_t> x(_num_points);
      std::vector<uint131_t> y(_num_points);

      HIP_CALL(hipMemcpy(x.data(), _dev_x, x.size() * sizeof(x[0]), hipMemcpyDeviceToHost));
      HIP_CALL(hipMemcpy(y.data(), _dev_y, y.size() * sizeof(y[0]), hipMemcpyDeviceToHost));

      for(int i = 0; i < 1000; i++) {
        uint131_t k1 = _priv_key_a[i];

        ecc::ecpoint_t sum = ecc::mul(k1, ecc::g());

        if(!ecc::exists(ecc::ecpoint_t(x[i], y[i]))) {
          LOG("{}", to_str(sum.x));
          LOG("{}", to_str(x[i]));
          LOG("Validation failed: Point does not exist");
          exit(1);
        }

        if(!ecc::is_equal(sum, ecc::ecpoint_t(x[i], y[i]))) {
          LOG("Expected:");
          LOG("{}", to_str(sum.x));
          LOG("{}", to_str(sum.y));
          LOG("{}", to_str(mont::from(sum.x)));
          LOG("{}", to_str(mont::from(sum.y)));
          
          LOG("Got:");
          LOG("{}", to_str(x[i]));
          LOG("{}", to_str(y[i]));
          LOG("{}", to_str(mont::from(x[i])));
          LOG("{}", to_str(mont::from(y[i])));

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
  HIP_IGNORE(hipFree(_staging_buf));
  HIP_IGNORE(hipFree(_staging_count));
  HIP_IGNORE(hipFree(_sanity_flag));

  HIP_IGNORE(hipFree(_walk_start));
}

void GPUPointFinder::allocate_buffers(int n)
{
  free_buffers();

  HIP_CALL(hipMalloc(&_mbuf, n * sizeof(uint131_t)));
  HIP_CALL(hipMalloc(&_dev_x, n * sizeof(uint131_t)));
  HIP_CALL(hipMalloc(&_dev_y, n * sizeof(uint131_t)));
  HIP_CALL(hipMalloc(&_dev_rx, _NUM_R_POINTS * sizeof(uint131_t)));
  HIP_CALL(hipMalloc(&_dev_ry, _NUM_R_POINTS * sizeof(uint131_t)));

  HIP_CALL(hipMallocManaged(&_priv_key_a, _num_points * sizeof(uint131_t)));
 
  HIP_CALL(hipMallocManaged(&_staging_buf, _staging_buf_size * sizeof(StagingPoint)));
  HIP_CALL(hipMallocManaged(&_staging_count, sizeof(uint32_t)));
  *_staging_count = 0;

  HIP_CALL(hipMallocManaged(&_result_buf, _result_buf_size * sizeof(DPResult)));

  HIP_CALL(hipMallocManaged(&_result_count, sizeof(uint32_t)));
  *_result_count = 0;
  
  HIP_CALL(hipMallocManaged(&_sanity_flag, sizeof(uint32_t)));
  *_sanity_flag = 0;


  HIP_CALL(hipMallocManaged(&_walk_start, _num_points * sizeof(uint64_t)));
}


void GPUPointFinder::step()
{

  // We don't need result buf for benchmark
  DPResult* result_buf = _benchmark ? nullptr : _result_buf;

  for(int i = 0; i < _iters_per_step; i++) {
    HIP_CALL(hipLaunchKernel(_do_step_ptr, dim3(_blocks), dim3(_threads), 0, _dev_x, _dev_y, _dev_rx, _dev_ry,
                  _mbuf, _num_points,
                  result_buf, _result_count,
                  _staging_buf, _staging_count,
                  _priv_key_a,
                  _counter,
                  _walk_start,
                  _dpmask));
    _counter++;
  }

  HIP_CALL(hipDeviceSynchronize());

  if(_first_run || _verify_points) {
    HIP_CALL(hipLaunchKernel((void*)_sanity_check_ptr, dim3(_blocks), dim3(_threads), 0, _dev_x, _dev_y, _num_points, _sanity_flag));
    HIP_CALL(hipDeviceSynchronize());

    if(*_sanity_flag != 0) {
      LOG("Detected {} errors\n", *_sanity_flag);

      std::vector<uint131_t> x(_num_points);
      std::vector<uint131_t> y(_num_points);

      HIP_CALL(hipMemcpy(x.data(), _dev_x, x.size() * sizeof(x[0]), hipMemcpyDeviceToHost));
      HIP_CALL(hipMemcpy(y.data(), _dev_y, y.size() * sizeof(y[0]), hipMemcpyDeviceToHost));

      int count = 0;
      for(int i = 0; i < _num_points; i++) {
        ecc::ecpoint_t p(x[i], y[i]);

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
    if(*_staging_count <= _staging_min) {
      refill_staging();
    }

    if(*_result_count >= _report_count) {
      report_points();
    }
  }
}

void GPUPointFinder::set_callback(std::function<void(const std::vector<DistinguishedPoint>&)> callback)
{
  _callback = callback;
}

size_t GPUPointFinder::work_per_step()
{
  return _num_points * 4;
}

int GPUPointFinder::iters_per_step()
{
  return _iters_per_step;
}


void GPUPointFinder::save_progress(const std::string& file_name)
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

  file.write((char*)&count, sizeof(count));
  if(!file.good()) {
    LOG("File not good");
    return;
  }

  // Write the counter
  file.write((char*)&_counter, sizeof(_counter));
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

  char* tmp = new char[sizeof(uint131_t) * count];

  // Points X
  HIP_CALL(hipMemcpy(tmp, _dev_x, sizeof(uint131_t) * count, hipMemcpyDeviceToHost));
  file.write(tmp, sizeof(uint131_t) * count);

  // Points Y
  HIP_CALL(hipMemcpy(tmp, _dev_y, sizeof(uint131_t) * count, hipMemcpyDeviceToHost));
  file.write(tmp, sizeof(uint131_t) * count);

  // Write starting counters
  memcpy(tmp, _walk_start, sizeof(uint64_t) * count);
  file.write(tmp, sizeof(uint64_t) * count);

  file.close();

  delete[] tmp;
}