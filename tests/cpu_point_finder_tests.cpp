#include <cassert>
#include <cstdio>
#include <string>
#include <vector>

#include "CPUPointFinder.h"
#include "CPUPointFinderF2N.h"
#include "ecc.h"

namespace {

void test_curve(const std::string& curve)
{
  ecc::set_curve(curve);

  auto finder = new CPUPointFinderF2N(1, 8);
  finder->init();
  finder->set_callback([](const std::vector<DistinguishedPoint>& points) {
    for(const DistinguishedPoint& point : points) {
      assert(ecc::exists(point.p));
      assert((point.p.x.w.v0 & 1) == 0);
    }
  });

  for(int i = 0; i < 16; i++)
    finder->step();

  const std::string progress = "/tmp/cpu_point_finder_test.dat";
  finder->save_progress(progress);

  auto loaded = new CPUPointFinderF2N(21);
  loaded->init(progress);
  assert(loaded->parallel_walks() == finder->parallel_walks());
  loaded->step();

  std::remove(progress.c_str());
}

} // namespace

int main()
{
  test_curve("ec2n79");
  return 0;
}
