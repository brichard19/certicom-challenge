#include <cassert>
#include <getopt.h>
#include <hip/hip_runtime.h>
#include "hip_helper.h"
#include <iostream>

#include "GPUPointFinder.h"
#include "ec_rho.h"
#include "util.h"

void benchmark(int hip_device)
{
  DistinguishedPointFinder* pf = new GPUPointFinder(hip_device, 0, 63, true);

  pf->init();

  util::Timer perf_timer;
  util::Timer run_timer;

  perf_timer.start();
  run_timer.start();
  size_t steps = 0;

  std::vector<double> ara;

  while(true) {
    pf->step();
    steps++;

    double t = perf_timer.elapsed();

    // Print performance info
    if(t >= 3.0) {
      size_t total = pf->work_per_step() * steps;

      double perf = (double)total / t;
      double iters = (double)steps * pf->iters_per_step() / t;

      perf_timer.start();

      steps = 0;

      std::cout << (perf / 1e6) << " MKeys/sec (" << iters << " iters/sec)" << std::endl;

      ara.push_back(perf);
    }

    // Exit after 30 seconds
    if(run_timer.elapsed() >= 15) {
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

int main(int argc, char**argv)
{
  int device = 0;

  while(true) {
    static struct option long_options[] = {
      {"gpu", required_argument, 0, 'g'},
    };

    int opt_idx = 0;

    int c = getopt_long(argc, argv, "g:", long_options, &opt_idx);

    if(c == -1) {
      break;
    }

    switch(c) {
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

  if(optind == argc) {
    std::cout << "Curve name required" << std::endl;
    return 1;
  }

  std::string curve_name = std::string(argv[optind]);

  if(curve_name != "ecp131" && curve_name != "ecp79") {
    std::cout << "Invalid curve name" << std::endl;
    return 1;
  }

  ecc::set_curve(curve_name);

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