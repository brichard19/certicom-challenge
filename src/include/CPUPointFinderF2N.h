#ifndef _CPU_POINT_FINDER_F2N_H
#define _CPU_POINT_FINDER_F2N_H

#include "CPUPointFinder.h"

class CPUPointFinderF2N : public CPUPointFinder {
private:
  bool _benchmark = false;
  std::vector<uint131_t> _chain;

public:
  explicit CPUPointFinderF2N(int dpbits, size_t num_points = 128, bool benchmark = false);
  void init() override;
  void init(const std::string& file) override;
  double step() override;

protected:
  void do_step_binary(std::vector<DistinguishedPoint>& results);
};

#endif
