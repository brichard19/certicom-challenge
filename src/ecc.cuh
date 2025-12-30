#ifndef _ECC_CUH
#define _ECC_CUH

#include <hip/hip_runtime.h>

#define POINTS_PER_THREAD 64

struct DPResult {
  uint131_t a;
  uint131_t x;
  uint131_t y;
  uint64_t length;
  uint32_t padding[12];
};

struct StagingPoint {
  uint131_t a;
  uint131_t x;
  uint131_t y;
};

extern "C" __global__ void clear_public_keys(uint131_t* x, uint131_t* y, int count);
extern "C" __global__ void sanity_check_p131(uint131_t* global_px, uint131_t* global_py, int count, int* errors);
extern "C" __global__ void sanity_check_p79(uint131_t* global_px, uint131_t* global_py, int count, int* errors);
extern "C" __global__ void sanity_check_p89(uint131_t* global_px, uint131_t* global_py, int count, int* errors);
extern "C" __global__ void reset_counters(uint64_t* start_pos, uint64_t value, int count);

extern "C" __global__ void do_step_p79(uint131_t* global_px, uint131_t* global_py, uint131_t* global_rx, uint131_t* global_ry, uint131_t* mbuf, int count, DPResult* result, int* result_count,
  StagingPoint* staging, int* staging_count, uint131_t* priv_key_a, uint64_t counter, uint64_t* start_pos, uint64_t dpmask);

  extern "C" __global__ void do_step_p89(uint131_t* global_px, uint131_t* global_py, uint131_t* global_rx, uint131_t* global_ry, uint131_t* mbuf, int count, DPResult* result, int* result_count,
  StagingPoint* staging, int* staging_count, uint131_t* priv_key_a, uint64_t counter, uint64_t* start_pos, uint64_t dpmask);

extern "C" __global__ void do_step_p131(uint131_t* global_px, uint131_t* global_py, uint131_t* global_rx, uint131_t* global_ry, uint131_t* mbuf, int count, DPResult* result, int* result_count,
  StagingPoint* staging, int* staging_count, uint131_t* priv_key_a, uint64_t counter, uint64_t* start_pos, uint64_t dpmask);

extern "C" __global__ void batch_multiply_p79(uint131_t* global_px, uint131_t* global_py, uint131_t* private_keys, uint131_t* mbuf, uint131_t* gx, uint131_t* gy, int priv_key_bit, int count);
extern "C" __global__ void batch_multiply_p89(uint131_t* global_px, uint131_t* global_py, uint131_t* private_keys, uint131_t* mbuf, uint131_t* gx, uint131_t* gy, int priv_key_bit, int count);
extern "C" __global__ void batch_multiply_p131(uint131_t* global_px, uint131_t* global_py, uint131_t* private_keys, uint131_t* mbuf, uint131_t* gx, uint131_t* gy, int priv_key_bit, int count);

#endif