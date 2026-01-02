#include "hip_helper.h"
#include <algorithm>
#include <cassert>
#include <getopt.h>
#include <hip/hip_runtime.h>
#include <iostream>

#include "GPUPointFinder.h"
#include "ec_rho.h"
#include "util.h"

double _benchmark_run_time = 10.0;

void benchmark(int hip_device)
{
  DistinguishedPointFinder* pf = new GPUPointFinder(hip_device, 63, true);

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

  delete pf;
}

int main(int argc, char** argv)
{
  std::string curve_name;
  int device = 0;

  while(true) {
    static struct option long_options[] = {
        {"curve", required_argument, 0, 'c'},
        {"gpu", required_argument, 0, 'g'},
    };

    int opt_idx = 0;

    int c = getopt_long(argc, argv, "c:g:", long_options, &opt_idx);

    if(c == -1) {
      break;
    }

    switch(c) {
    case 'c':
      curve_name = std::string(optarg);
      break;

    case 'g':
      device = atoi(optarg);
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

  // Check device ID
  int device_count = 0;
  HIP_CALL(hipGetDeviceCount(&device_count));

  if(device >= device_count) {
    std::cout << "Invalid device " << device << std::endl;
    return 1;
  }

  benchmark(device);

  return 0;
}