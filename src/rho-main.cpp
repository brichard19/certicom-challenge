#if defined(BUILD_GPU) && defined(BUILD_CPU)
#error "BUILD_GPU and BUILD_CPU cannot both be defined"
#elif !defined(BUILD_GPU) && !defined(BUILD_CPU)
#error "One of BUILD_GPU or BUILD_CPU must be defined"
#endif

#include "fmt/format.h"
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

#include "binary_encoder.h"
#include "ec_rho.h"
#include "log.h"
#include "montgomery.h"
#include "signal_handler.h"
#include "util.h"

#if defined(BUILD_GPU)
#include "GPUPointFinder.h"
#include "hip_helper.h"
#include <hip/hip_runtime.h>
#endif

#if defined(BUILD_CPU)
#include "CPUPointFinder.h"
#include "CPUPointFinderF2N.h"
#endif

#if defined(BUILD_MPI)
#include "mpi_helper.h"
#include <mpi.h>
#endif

namespace {

const double _save_interval = 60.0;
const double _perf_interval = 5.0;

volatile bool _running = true;

std::string _data_file = "";
std::string _data_dir = "";
std::string _results_dir = "";
std::string _hostname = "";

// MPI variables
bool _use_mpi = false;
int _world_size = -1;
int _world_rank = -1;
volatile bool _mpi_thread_running = true;

int _dpbits = 0;
std::string _curve_name;

#if defined(BUILD_GPU)
HIPDeviceMap _device_map;
#endif

// Used for both CPU and GPU
int _device_idx = 0;

} // namespace

std::string get_data_file_name()
{
#if defined(BUILD_GPU)
  return fmt::format("GPU-{}.dat", _device_map[_device_idx].uuid);
#elif defined(BUILD_CPU)
  return fmt::format("CPU-{}.dat", _device_idx);
#endif
}

DistinguishedPointFinder* create_point_finder()
{
#if defined(BUILD_GPU)
  return new GPUPointFinder(_device_idx, _dpbits);
#elif defined(BUILD_CPU)
  return new CPUPointFinderF2N(_dpbits, 256);
#endif
}

// Saves distingusihed points to disk
void save_to_disk(const std::vector<DistinguishedPoint>& dps)
{
  std::random_device random;
  std::uniform_int_distribution<uint64_t> distribution;
  uint64_t high = distribution(random);
  uint64_t low = distribution(random);
  std::string id = fmt::format("{:016x}{:016x}", high, low);

  // Use a temp file when writing since another thread periodically picks up all the
  // .dat files
  std::string tmp_name = fmt::format("{}/{}.tmp", _results_dir, id);
  std::string file_name = fmt::format("{}/{}.dat", _results_dir, id);

  auto encoded = encode_dps(dps, ecc::curve_strength(), _dpbits);
  std::ofstream of(tmp_name, std::ios::binary);

  of.write((const char*)encoded.data(), encoded.size());
  of.close();

  std::filesystem::rename(tmp_name, file_name);
}

void dp_callback(const std::vector<DistinguishedPoint>& dps)
{
  uint64_t dpmask = ((uint64_t)1 << _dpbits) - 1;

  LOG("Found {} distinguished points:", dps.size());

  // Validate they are correct
  for(auto dp : dps) {
    assert(ecc::exists(dp.p));
    assert((dp.p.x.w.v0 & dpmask) == 0);
  }

  save_to_disk(dps);
}

void main_loop()
{

  std::string data_file_path = _data_dir + "/" + _hostname + "/" + _data_file;

  DistinguishedPointFinder* pf = create_point_finder();

  pf->init(data_file_path);

  pf->set_callback(dp_callback);

  util::Timer perf_timer;
  util::Timer save_timer;

  perf_timer.start();
  save_timer.start();
  size_t steps = 0;
  double accu_time = 0.0;

  while(_running) {
    accu_time += pf->step();
    steps++;

    double t = perf_timer.elapsed();

    // Print performance info
    if(t >= _perf_interval) {
      size_t total = pf->work_per_step() * steps;

      double perf = (double)total / accu_time;
      double iters = (double)steps * pf->iters_per_step() / accu_time;

      perf_timer.start();

      LOG("Perf: {:.2f} MKeys/sec | Iters: {:.2f} iters/sec | Parallel walks: {}", perf / 1e6,
          iters, pf->parallel_walks());

      steps = 0;
      accu_time = 0;
    }

    // Save data
    t = save_timer.elapsed();
    if(t >= _save_interval) {
      pf->save_progress(data_file_path);
      save_timer.start();
    }
  }

  // Save data
  pf->save_progress(data_file_path);

  delete pf;
}

void signal_handler(int signal)
{
  std::cout << "Exiting..." << std::endl;
  _running = false;
}

bool init_directories()
{
  _results_dir = _data_dir + "/results";

  std::error_code err;

  if(!std::filesystem::create_directories(_data_dir, err) && !std::filesystem::exists(_data_dir)) {
    std::cout << fmt::format("Error creating directory '{}': {}", _data_dir, err.message())
              << std::endl;
    return false;
  }

  if(!std::filesystem::create_directories(_results_dir, err) &&
     !std::filesystem::exists(_results_dir)) {
    std::cout << fmt::format("Error creating directory '{}': {}", _results_dir, err.message())
              << std::endl;
    return false;
  }

  std::string progress_dir = _data_dir + "/" + _hostname;
  if(!std::filesystem::create_directories(progress_dir, err) &&
     !std::filesystem::exists(progress_dir)) {
    std::cout << fmt::format("Error creating directory '{}': {}", progress_dir, err.message())
              << std::endl;
    return false;
  }

  return true;
}

int main(int argc, char** argv)
{

  _hostname = util::get_hostname();

  bool gpu_flag = false;

  while(true) {
    static struct option long_options[] = {
#if defined(BUILD_GPU)
        {"gpu", required_argument, 0, 'g'},
#endif
        {"data-dir", required_argument, 0, 'd'},
        {"mpi", no_argument, 0, 'm'},
        {"curve", required_argument, 0, 'c'},
        {"dp-bits", required_argument, 0, 'b'},
        {NULL, 0, NULL, 0}};

    int opt_idx = 0;

    int c = getopt_long(argc, argv, "", long_options, &opt_idx);

    if(c == -1) {
      break;
    }

    switch(c) {
    case 'c':
      _curve_name = std::string(optarg);
      break;

    case 'd':
      _data_dir = std::string(optarg);
      break;
#if defined(BUILD_GPU)
    case 'g':
      _device_idx = std::stoi(optarg);
      gpu_flag = true;
      break;
#endif
    case 'b':
      _dpbits = std::stoi(optarg);
      break;

    case 'm':
#if defined(BUILD_MPI)
      _use_mpi = true;
#else
      std::cout << "Error: using --mpi but not built with MPI support!" << std::endl;
      return 1;
#endif
      break;

    case '?':
      return 1;
      break;

    default:
      std::cout << "Invalid argument" << std::endl;
      return 1;
    }
  }

#if defined(BUILD_GPU)
  _device_map = get_device_map();
#endif

  if(_curve_name.empty()) {
    std::cout << "--curve required" << std::endl;
    return 1;
  }

  if(_data_dir.empty()) {
    std::cout << "--data required" << std::endl;
    return 1;
  }

  if(_dpbits == 0) {
    std::cout << "--dp-bits required" << std::endl;
    return 1;
  }

  if(_dpbits < 16 || _dpbits > 31) {
    std::cout << "--dp-bits must be between 16 and 30" << std::endl;
    return 1;
  }

  try {
    ecc::set_curve(_curve_name);
#if defined(BUILD_GPU)
    if(ecc::is_prime_curve() == false) {
      std::cout << "Curve " << _curve_name << " is not supported on GPU" << std::endl;
      return 1;
    }
#endif
#if defined(BUILD_CPU)
    if(ecc::is_binary_curve() == false) {
      std::cout << "Curve " << _curve_name << " is not supported on CPU" << std::endl;
      return 1;
    }
#endif
  } catch(...) {
    std::cout << "Invalid curve name" << std::endl;
    return 1;
  }

#if defined(BUILD_GPU)
  if(ecc::is_prime_curve() == false) {
    std::cout << "GPU builds only support prime curves" << std::endl;
    return 1;
  }
#elif defined(BUILD_CPU)
  if(ecc::is_binary_curve() == false) {
    std::cout << "CPU builds only support binary curves" << std::endl;
    return 1;
  }
#endif

  if(_use_mpi && gpu_flag) {
    std::cout << "-g flag incompable when using MPI" << std::endl;
    return 1;
  }

#if defined(BUILD_GPU)
  // Check device ID
  int device_count = 0;
  HIP_CALL(hipGetDeviceCount(&device_count));

  if(device_count == 0) {
    std::cout << "No GPUs available" << std::endl;
    return 1;
  }
#endif

#if defined(BUILD_MPI)

  // Initialize MPI, select device
  if(_use_mpi) {
    if(MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      std::cout << " MPI_Init failed" << std::endl;
      return 1;
    }

    MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &_world_size));

    MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &_world_rank));

    // TODO: Is this the correct way to get the local rank?
    int local_rank = -1;
    MPI_Comm local_comm;
    MPI_CALL(MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, _world_rank, MPI_INFO_NULL,
                                 &local_comm));
    MPI_CALL(MPI_Comm_rank(local_comm, &local_rank));
    MPI_CALL(MPI_Comm_free(&local_comm));

#if defined(BUILD_GPU)
    // One GPU per rank
    if(local_rank >= device_count) {
      MPI_Finalize();
      return 1;
    }
#endif

    _device_idx = local_rank;
  }
#endif

#if defined(BUILD_GPU)
  if(_device_idx >= device_count) {
    std::cout << "Invalid device " << _device_idx << std::endl;
    return 1;
  }
#endif

  if(!init_directories()) {
    return 1;
  }

  // Use device UUID for data file
  _data_file = get_data_file_name();

  // Set interrupt handler
  set_signal_handler(signal_handler);

  // Run main loop
  main_loop();

#if defined(BUILD_MPI)
  // Cleanup MPI
  if(_use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
