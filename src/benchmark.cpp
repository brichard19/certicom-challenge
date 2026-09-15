#include <algorithm>
#include <cassert>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <optional>

#if defined(BUILD_GPU)
#include "GPUPointFinder.h"
#include "hip_helper.h"
#include <hip/hip_runtime.h>
#endif

#if defined(BUILD_CPU)
#include "CPUPointFinder.h"
#include "CPUPointFinderF2N.h"
#endif

#include "ec_rho.h"
#include "util.h"

double _benchmark_run_time = 10.0;

int _hip_device = 0;

DistinguishedPointFinder* create_point_finder()
{
#if defined(BUILD_GPU)
  return new GPUPointFinder(_hip_device, 20, true);
#elif defined(BUILD_CPU)
  return new CPUPointFinderF2N(20, 256, true);
#else
#error "Either BUILD_GPU or BUILD_CPU must be defined"
#endif
}

void benchmark()
{
  DistinguishedPointFinder* pf = create_point_finder();

  pf->init();

  util::Timer run_timer;

  run_timer.start();
  size_t steps = 0;

  std::vector<double> ara;

  double gpu_time = 0.0;

  double total_time = 0.0;

  while(true) {
    gpu_time += pf->step();
    steps++;

    // Print performance info
    double t = run_timer.elapsed();
    if(t >= 3.0) {
      total_time += t;
      run_timer.start();

      size_t total = pf->work_per_step() * steps;
      double perf = (double)total / gpu_time;
      double iters = (double)steps * pf->iters_per_step() / gpu_time;

      steps = 0;
      gpu_time = 0;
      std::cout << (perf / 1e6) << " MKeys/sec (" << iters << " iters/sec)" << std::endl;

      ara.push_back(perf);
    }

    if(total_time >= _benchmark_run_time) {
      break;
    }
  }

  delete pf;

  // Remove lowest value then take the average
  std::sort(ara.begin(), ara.end());

  // Take average. Ignore the lowest.
  double sum = 0.0;
  for(int i = 1; i < ara.size(); i++) {
    sum += ara[i];
  }

  double avg = sum / (ara.size() - 1);

  std::cout << std::endl;
  std::cout << (avg / 1e6) << " MKeys/sec" << std::endl;
}

int main(int argc, char** argv)
{
  std::string curve_name;

  while(true) {
    static struct option long_options[] = {
#if defined(BUILD_GPU)
        {"gpu", required_argument, 0, 'g'},
#endif
        {"curve", required_argument, 0, 'c'},
        {"cpu", no_argument, 0, 'C'},
        {NULL, 0, NULL, 0},
    };

    int opt_idx = 0;

    int c = getopt_long(argc, argv, "c:g:C", long_options, &opt_idx);

    if(c == -1) {
      break;
    }

    switch(c) {
    case 'c':
      curve_name = std::string(optarg);
      break;

    case 'g':
      _hip_device = atoi(optarg);
      break;

    case 'C':
      // No action needed, as the default is CPU
      break;

    case '?':
      break;

    default:
      std::cout << "Invalid argument" << std::endl;
      exit(1);
    }
  }

  if(curve_name.empty()) {
    std::cout << "--curve required" << std::endl;
    return 1;
  }

  try {
    ecc::set_curve(curve_name);
  } catch(...) {
    std::cout << "Invalid curve name" << std::endl;
    return 1;
  }

#if defined(BUILD_GPU)
  int device_count = 0;
  HIP_CALL(hipGetDeviceCount(&device_count));

  if(_hip_device < 0 || _hip_device >= device_count) {
    std::cout << "Invalid device " << _hip_device << std::endl;
    return 1;
  }
#endif

  benchmark();

  return 0;
}
