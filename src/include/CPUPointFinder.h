#ifndef _CPU_POINT_FINDER_H
#define _CPU_POINT_FINDER_H

#include <functional>
#include <string>
#include <vector>

#include "ec_rho.h"

class CPUPointFinder : public DistinguishedPointFinder {
private:
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

  std::vector<uint131_t> _chain;

  std::function<void(const std::vector<DistinguishedPoint>&)> _callback;

  void generate_walk(size_t index);
  void step_prime(std::vector<DistinguishedPoint>& results);
  void step_binary(std::vector<DistinguishedPoint>& results);
  void load(const std::string& file);

public:
  explicit CPUPointFinder(int dpbits, size_t num_points = 128);

  void init() override;
  void init(const std::string& file) override;
  double step() override;
  void set_callback(std::function<void(const std::vector<DistinguishedPoint>&)> callback) override;
  void save_progress(const std::string& file) override;
  size_t work_per_step() override;
  int iters_per_step() override;
  int parallel_walks() override;
};

#endif
