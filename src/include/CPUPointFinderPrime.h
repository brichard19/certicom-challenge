#ifndef _CPU_POINT_FINDER_PRIME_H
#define _CPU_POINT_FINDER_PRIME_H

#include "CPUPointFinder.h"

class CPUPointFinderPrime : public CPUPointFinder {
private:
  std::vector<uint131_t> _chain;

public:
  explicit CPUPointFinderPrime(int dpbits, size_t num_points = 128);
  void init() override;
  void init(const std::string& file) override;
  double step() override;

protected:
  void do_step_prime(std::vector<DistinguishedPoint>& results);
};

#endif
