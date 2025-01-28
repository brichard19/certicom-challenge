#ifndef _GPU_POINT_FINDER_H
#define _GPU_POINT_FINDER_H

#include <vector>

#include "ec_rho.h"
#include "ecc.cuh"

class GPUPointFinder : public DistinguishedPointFinder {

private:
#ifdef DEBUG
  const bool _verify_points = true;
#else
  const bool _verify_points = false;
#endif

  void* _do_step_ptr = nullptr;
  void* _batch_multiply_ptr = nullptr;
  void* _sanity_check_ptr = nullptr;

  bool _first_run = true;

  uint32_t _num_points = 1024 * 1024;
  int _NUM_R_POINTS = 32;

  int _iters_per_step = 4;
  int _result_buf_size;
  int _staging_buf_size = 1024;
  int _report_count = 16;
  int _staging_min = 16;

  int _dpbits = 16;
  uint64_t _dpmask;

  std::function<void(const std::vector<DistinguishedPoint>&)> _callback;

  int _device;
  uint131_t* _mbuf = nullptr;
  uint131_t* _dev_x = nullptr;
  uint131_t* _dev_y = nullptr;
  uint131_t* _dev_rx = nullptr;
  uint131_t* _dev_ry = nullptr;

  uint131_t* _priv_key_a = nullptr;

  StagingPoint* _staging_buf = nullptr;
  uint32_t* _staging_count = nullptr;

  DPResult* _result_buf = nullptr;
  uint32_t* _result_count = nullptr;

  uint32_t* _sanity_flag = nullptr;

  uint64_t _counter = 0;
  uint64_t* _walk_start = nullptr;

  int _blocks = 1;
  int _threads = 32;

  std::vector<uint131_t> _rx;
  std::vector<uint131_t> _ry;

  void load(const std::string& file_name);

  void allocate_buffers(int n);

  void free_buffers();

  void refill_staging();

  void report_points();

public:

  GPUPointFinder(int device, uint32_t num_points, int dpbits);

  virtual ~GPUPointFinder();

  virtual void init();

  virtual void init(const std::string& file);

  virtual void step();

  virtual void set_callback(std::function<void(const std::vector<DistinguishedPoint>&)> callback);

  virtual void save_progress(const std::string& file);

  virtual size_t work_per_step();
  
  virtual int iters_per_step();
};

#endif