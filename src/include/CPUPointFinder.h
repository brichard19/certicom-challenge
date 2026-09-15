#ifndef _CPU_POINT_FINDER_H
#define _CPU_POINT_FINDER_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "ec_rho.h"

class CPUPointFinder : public DistinguishedPointFinder {
protected:
  int _dpbits;
  uint64_t _dpmask;
  size_t _num_points;
  int _iters_per_step = 1;
  uint64_t _counter = 0;

  std::vector<uint131_t> _rx;
  std::vector<uint131_t> _ry;

  std::vector<uint131_t> _x;
  std::vector<uint131_t> _y;
  std::vector<uint131_t> _priv;
  std::vector<uint64_t> _walk_len;

  std::function<void(const std::vector<DistinguishedPoint>&)> _callback;

  void generate_walk(size_t index);
  void load(const std::string& file);
  void save_progress_impl(const std::string& file);

public:
  explicit CPUPointFinder(int dpbits, size_t num_points = 128);

  // Abstract interface: derived classes must implement these
  void init() override = 0;
  void init(const std::string& file) override = 0;
  double step() override = 0;

  virtual ~CPUPointFinder();

  // Common implementations
  void set_callback(std::function<void(const std::vector<DistinguishedPoint>&)> callback) override;
  void save_progress(const std::string& file) override;
  size_t work_per_step() override;
  int iters_per_step() override;
  int parallel_walks() override;
};

#endif
