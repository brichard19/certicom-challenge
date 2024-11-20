#include <cassert>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <hip/hip_runtime.h>
#include <iostream>
#include <mutex>
#include <thread>
#include <unistd.h>

#include "GPUPointFinder.h"
#include "binary_encoder.h"
#include "ec_rho.h"
#include "fmt/format.h"
#include "hip_helper.h"
#include "http_client.h"
#include "log.h"
#include "signal_handler.h"
#include "util.h"


#ifdef USE_MPI

#include <mpi.h>

#define MPI_CALL(condition)\
{\
  int err = condition;\
  if(err != MPI_SUCCESS) {\
    std::cout << "MPI error " << err << " " << __LINE__ << std::endl;\
    MPI_Finalize();\
    exit(1);\
  }\
}\

#endif


#define DEFAULT_FILE_FORMAT_STRING "{}-gpu{}.dat"
#define DEFAULT_DATA_DIRECTORY "data"

namespace {

  const double _save_interval = 600.0;
  const double _perf_interval = 5.0;

#if defined(CURVE_P131)
  const int _dpbits = 32;
#elif defined(CURVE_P79)
  const int _dpbits = 28;
#else
#error "Curve is undefined"
#endif

  std::string _url = "127.0.0.1";
  int _port = 8080;

  volatile bool _running = true;

  int _hip_device = 0;

  std::string _data_file = "";
  std::string _data_dir = DEFAULT_DATA_DIRECTORY;
  std::string _results_dir;

  // System will choose number of points 
  uint32_t _num_points = 0;

  // MPI variables
  bool _use_mpi = false;
  int _world_size;
  int _world_rank;
}

// Saves distingusihed points to disk
void save_to_disk(const std::vector<DistinguishedPoint>& dps)
{
  // Use a temp file when writing since another thread periodically picks up all the
  // .txt files
  std::string tmp_name = fmt::format("{}/{}.tmp", _results_dir, (int)time(NULL));
  std::string file_name = fmt::format("{}/{}.txt", _results_dir, (int)time(NULL));

  std::ofstream of(tmp_name);

  for(auto p : dps) {
    of << fmt::format("{} {} {} {}", to_str(p.a), to_str(p.p.x), to_str(p.p.y), p.length) << std::endl;
  }
  of.close();

  std::filesystem::rename(tmp_name, file_name);
}

void callback(const std::vector<DistinguishedPoint>& dps)
{
  uint64_t dpmask = ((uint64_t)1 << _dpbits) - 1;

  LOG("Found {} distinguished points:", dps.size());

  for(auto dp : dps) {
    assert(ecc::exists(dp.p));
    assert((dp.p.x.v[0] & dpmask) == 0);
    //LOG("{} {} {}", to_str(dp.p.x), to_str(dp.p.y), dp.length);
  }

  save_to_disk(dps);
}



// Encodes a DistinguishedPoint into a string of bytes
std::vector<uint8_t> encode_dp(const DistinguishedPoint& dp)
{
  // 17 + 17 + 8 bytes
  std::vector<uint8_t> buf(42);
  uint8_t* ptr = buf.data();

  // x-coordinate, 17 bytes
  // TODO: Compress further by removing the "distinguished bits", which
  // are all 0's. 
  memcpy(ptr, dp.p.x.v, 17);
  
  // Append sign bit to x coordinate. The x coordinate is 131 bits,
  // so set the 132nd bit.
  uint8_t sign = is_odd(dp.p.y) ? 1 : 0;
  ptr[16] |= (sign << 3);

  ptr += 17;

  // a-exponent, 17 bytes
  memcpy(ptr, dp.a.v, 17);
  ptr += 17;

  // Walk length, 8 bytes
  memcpy(ptr, &dp.length, sizeof(uint64_t));

  // TODO: User ID?

  return buf;
}

void results_thread_function()
{
  bool defer_upload = false;

  util::Timer defer_timer;

  while(_running) {

    sleep(5);

    if(defer_upload && defer_timer.elapsed() < 600.0) {
      continue;
    } else {
      defer_upload = false;
    }

    std::vector<std::string> files;

    // Get list of results files
    for(const auto& entry : std::filesystem::directory_iterator(_results_dir)) {
      std::string f = entry.path().filename();
      if(f.rfind(".txt") != std::string::npos) {
        files.push_back(entry.path().filename());
      }
    }

    if(files.size() == 0) {
      continue;
    }

    // Load results from files
    std::vector<DistinguishedPoint> points;

    for(auto& file : files) {
      std::ifstream f(_results_dir + "/" + file, std::ios::in);

      std::string a_hex;
      std::string x_hex;
      std::string y_hex;
      uint64_t length;

      while(f >> a_hex >> x_hex >> y_hex >> length) {
        uint131_t a = make_uint131(a_hex);
        uint131_t x = make_uint131(x_hex);
        uint131_t y = make_uint131(y_hex);

        DistinguishedPoint dp(a, ecc::ecpoint_t(x, y), length);
        points.push_back(dp);
      }
    }

    // Encode results 
    bool success = true;

    BinaryEncoder encoder;

    uint32_t num_points = points.size();

    encoder.encode(num_points);

    // TODO: Encode points
    for(auto p : points) {
      std::vector<uint8_t> buf = encode_dp(p);
      encoder.encode(buf.data(), buf.size());
    }
    
    LOG("Uploading {} points to server", points.size());

    // TODO: Use base64 instead (smaller data)
    // Upload to server 
    std::string hex = util::to_hex(encoder.get_ptr(), encoder.get_size());
    //LOG("{}", hex);

    try {
      HTTPClient http(_url, _port);
      http.submit_points(hex);
      defer_upload = false;
    }catch(std::exception& ex) {
      LOG("Upload error: {}", ex.what());
      success = false;
      
      LOG("Deferring upload for 10 minutes");
      defer_upload = true;
      defer_timer.start();
    }

    // Delete files after successful upload
    if(success) {
      for(auto& f : files) {
        std::filesystem::remove(_results_dir + "/" + f);
      }
    }
  }
}

void main_loop()
{

  std::string data_file_path = _data_dir + "/" + _data_file;

  DistinguishedPointFinder* pf = new GPUPointFinder(_hip_device, _num_points, _dpbits);

  pf->init(data_file_path);

  pf->set_callback(callback);

  util::Timer perf_timer;
  util::Timer save_timer;

  perf_timer.start();
  save_timer.start();
  size_t steps = 0;

  while(_running) {
    pf->step();
    steps++;

    double t = perf_timer.elapsed();

    // Print performance info
    if(t >= _perf_interval) {
      size_t total = pf->work_per_step() * steps;

      double perf = (double)total / t;
      double iters = (double)steps / t;

      perf_timer.start();

      steps = 0;


      LOG("{:.2f} MKeys/sec ({:.2f} iters/sec)", perf / 1e6, iters);
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

  if(!std::filesystem::create_directories(_data_dir, err)) {
    if(std::filesystem::exists(_data_dir)) {
      return true;
    }
    std::cout << fmt::format("Error creating directory '{}': {}", _data_dir, err.message()) << std::endl;
    return false;
  }

  if(!std::filesystem::create_directories(_results_dir, err)) {
    if(std::filesystem::exists(_data_dir)) {
      return true;
    }
    std::cout << fmt::format("Error creating directory '{}': {}", _data_dir, err.message()) << std::endl;
    return false;
  }

  return true;
}

int main(int argc, char**argv)
{

#if defined(CURVE_P131)
  ecc::set_curve("ecp131");
#elif defined(CURVE_P79)
  ecc::set_curve("ecp79");
#else
#error "Curve is undefined"
#endif

  bool gpu_flag = false;

  while(true) {
    static struct option long_options[] = {
      {"gpu", required_argument, 0, 'g'},
      {"data", required_argument, 0, 'd'},
      {"file", required_argument, 0, 'f'},
      {"mpi", no_argument, 0, 'm'},
      {"size", required_argument, 0, 's'},
      {NULL, 0, NULL, 0}
    };

    int opt_idx = 0;

    int c = getopt_long(argc, argv, "d:f:g:s:m", long_options, &opt_idx);

    if(c == -1) {
      break;
    }

    switch(c) {
      case 'd':
        _data_dir = std::string(optarg);
        break;

      case 'f':
        _data_file = std::string(optarg);
        break;

      case 'g':
        _hip_device = atoi(optarg);
        gpu_flag = true;
        break;

      case 'm':
        _use_mpi = true;
        break;

      case 's':
        _num_points = atoi(optarg);
        break;

      case '?':
        break;

      default:
        std::cout << "Invalid argument" << std::endl;
        exit(1);
    }
  }

  if(_use_mpi && gpu_flag) {
    std::cout << "-g flag incompable when using MPI" << std::endl;
    return 1;
  }

  // Check device ID
  int device_count = 0;
  HIP_CALL(hipGetDeviceCount(&device_count));

  if(device_count == 0) {
    std::cout << "No GPUs available" << std::endl;
    return 1;
  }


#ifdef USE_MPI

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
    MPI_CALL(MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, _world_rank, MPI_INFO_NULL, &local_comm));
    MPI_CALL(MPI_Comm_rank(local_comm, &local_rank));
    MPI_CALL(MPI_Comm_free(&local_comm));

    // One GPU per rank
    if(local_rank >= device_count) {
      MPI_Finalize();
      return 1;
    }
 
    _hip_device = local_rank;
  }
#endif

  if(_hip_device >= device_count) {
    std::cout << "Invalid device " << _hip_device << std::endl;
    return 1;
  }

  if(!init_directories()) {
    return 1;
  }

  // Use default file unless one was provided
  if(_data_file.empty()) {
    _data_file = fmt::format(DEFAULT_FILE_FORMAT_STRING, ecc::curve_name().c_str(), _hip_device);
  }

  // Set interrupt handler
  set_signal_handler(signal_handler);

  // Start point reporter thread
  std::thread results_thread(results_thread_function);

  // Run main loop
  main_loop();

  // Wait for point thread to finish
  results_thread.join();

#ifdef USE_MPI
  // Cleanup MPI
  if(_use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}