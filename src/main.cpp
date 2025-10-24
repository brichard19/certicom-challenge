#include <cassert>
#include <chrono>
#include <filesystem>
#include "fmt/format.h"
#include <fstream>
#include <getopt.h>
#include <hip/hip_runtime.h>
#include <iostream>
#include <mutex>
#include <thread>
#include <map>
#include <vector>

#include "GPUPointFinder.h"
#include "binary_encoder.h"
#include "ec_rho.h"
#include "hip_helper.h"
#include "log.h"
#include "signal_handler.h"
#include "util.h"
#include "montgomery.h"

#ifdef BUILD_MPI
#include <mpi.h>
#include "mpi_helper.h"
#endif

#ifdef BUILD_HTTP
#include "http_client.h"
#endif

namespace {

  std::map<std::string, int> _curve_bits = {
    {"ecp131", 131},
    {"ecp79", 79}
  };

  const double _save_interval = 600.0;
  const double _perf_interval = 5.0;

  int _port = 8080;

  volatile bool _running = true;

  HIPDeviceMap _device_map;
  int _hip_device = 0;

  std::string _data_file = "";
  std::string _data_dir = "";
  std::string _results_dir = "";
  std::string _hostname = "";

  // System will choose number of points 
  uint32_t _num_points = 0;

  // MPI variables
  bool _use_mpi = false;
  int _world_size = -1;
  int _world_rank = -1;
  volatile bool _mpi_thread_running = true;

  int _dpbits = 32;
  std::string _curve_name;
  std::string _upload_server = "";
  bool _use_upload = false;
}

// Saves distingusihed points to disk
void save_to_disk(const std::vector<DistinguishedPoint>& dps)
{
  // Use a temp file when writing since another thread periodically picks up all the
  // .dat files
  //std::string tmp_name = fmt::format("{}/{}.tmp", _results_dir, (int)time(NULL));
  std::string tmp_name = fmt::format("{}/{}.tmp", _results_dir, (int)time(NULL));
  std::string file_name = fmt::format("{}/{}.dat", _results_dir, (int)time(NULL));

  auto encoded = encode_dps(dps, _dpbits, _curve_bits[_curve_name]);
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

  if(_use_mpi == false || (_use_mpi == true && _world_rank == 0)) {
    save_to_disk(dps);
  } else {
#ifdef BUILD_MPI
    LOG("Rank {} reporting {} points", _world_rank, dps.size());
    MPI_CALL(MPI_Send(dps.data(), dps.size() * sizeof(dps[0]), MPI_BYTE, 0, 0, MPI_COMM_WORLD));
#endif

  }
}



#ifdef BUILD_MPI
void mpi_recv_thread_function()
{
  int buf_size = 128 * 1024;
  std::vector<char> buf(buf_size);
  
  LOG("MPI thread started");

  MPI_Request request;

  // Async request
  MPI_CALL(MPI_Irecv(buf.data(), buf_size, MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &request));
  
  while(_mpi_thread_running == true) {
    int flag = 0;
    MPI_Status status;
    MPI_CALL(MPI_Test(&request, &flag, &status));
    if(flag) {
      int num_bytes;
      MPI_CALL(MPI_Get_count(&status, MPI_BYTE, &num_bytes));

      printf("MPI: Received %d bytes\n", num_bytes);

      assert(num_bytes % sizeof(DistinguishedPoint) == 0);

      int num_points = num_bytes / sizeof(DistinguishedPoint);
      std::vector<DistinguishedPoint> dps(num_points);

      memcpy(dps.data(), buf.data(), num_bytes);

      save_to_disk(dps);
  
      // New async request
      MPI_CALL(MPI_Irecv(buf.data(), buf_size, MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &request));
    }
    std::this_thread::sleep_for(std::chrono::seconds(3));
  }
}
#endif

#ifdef BUILD_HTTP
void upload_thread_function()
{
  bool defer_upload = false;

  util::Timer defer_timer;

  while(_running) {

    std::this_thread::sleep_for(std::chrono::seconds(5));

    if(_running && defer_upload && defer_timer.elapsed() < 600.0) {
      continue;
    } else {
      defer_upload = false;
    }

    std::vector<std::string> files;

    // Get list of results files
    for(const auto& entry : std::filesystem::directory_iterator(_results_dir)) {
      std::string f = entry.path().filename().string();
      if(f.rfind(".dat") != std::string::npos) {
        files.push_back(f);
      }
    }

    if(files.size() == 0) {
      continue;
    }

    // Load results from files
    std::vector<DistinguishedPoint> points;

    for(auto& file : files) {
      std::ifstream f(_results_dir + "/" + file, std::ios::binary);

      f.seekg(0, std::ios::end);
      std::streamsize file_size = f.tellg();
      f.seekg(0, std::ios::beg);

      std::vector<uint8_t> buf(file_size);
      f.read((char*)buf.data(), buf.size());

      auto dps = decode_dps(buf.data(), buf.size());

      points.insert(points.end(), dps.begin(), dps.end());      
    }

    // Encode results 
    bool success = true;

    auto encoded = encode_dps(points, _dpbits, _curve_bits[_curve_name]);

    LOG("Uploading {} points to server", points.size());

    // TODO: Use base64 instead (smaller data)
    // Upload to server 
    std::string hex = util::to_hex(encoded.data(), encoded.size());

    try {
      HTTPClient http(_upload_server, _port);
      http.submit_points(_curve_name, hex);
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

#endif

void main_loop()
{

  std::string data_file_path = _data_dir + "/" + _hostname + "/" + _data_file;

  DistinguishedPointFinder* pf = new GPUPointFinder(_hip_device, _num_points, _dpbits);

  pf->init(data_file_path);

  pf->set_callback(dp_callback);

  util::Timer perf_timer;
  util::Timer save_timer;

  perf_timer.start();
  save_timer.start();
  size_t steps = 0;
  double gpu_time = 0.0;

  while(_running) {
    gpu_time += pf->step();
    steps++;

    double t = perf_timer.elapsed();

    // Print performance info
    if(t >= _perf_interval) {
      size_t total = pf->work_per_step() * steps;

      double perf = (double)total / gpu_time;
      double iters = (double)steps * pf->iters_per_step() / gpu_time;

      perf_timer.start();

      LOG("{:.2f} MKeys/sec ({:.2f} iters/sec)", perf / 1e6, iters);

      steps = 0;
      gpu_time = 0;
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
    std::cout << fmt::format("Error creating directory '{}': {}", _data_dir, err.message()) << std::endl;
    return false;
  }

  if(!std::filesystem::create_directories(_results_dir, err) && !std::filesystem::exists(_results_dir)) {
    std::cout << fmt::format("Error creating directory '{}': {}", _results_dir, err.message()) << std::endl;
    return false;
  }

  std::string progress_dir = _data_dir + "/" + _hostname;
  if(!std::filesystem::create_directories(progress_dir, err) && !std::filesystem::exists(progress_dir)) {
    std::cout << fmt::format("Error creating directory '{}': {}", progress_dir, err.message()) << std::endl;
    return false;
  }

  return true;
}

int main(int argc, char**argv)
{
 
  _hostname = util::get_hostname();

  bool gpu_flag = false;
 
  while(true) {
    static struct option long_options[] = {
      {"gpu", required_argument, 0, 'g'},
      {"data-dir", required_argument, 0, 'd'},
      {"file", required_argument, 0, 'f'},
      {"mpi", no_argument, 0, 'm'},
      {"size", required_argument, 0, 's'},
      {"curve", required_argument, 0, 'c'},
      {"upload-server", required_argument, 0, 'u'},
      {NULL, 0, NULL, 0}
    };

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
      
      case 'u':
        _use_upload = true;
        _upload_server = std::string(optarg);
        break;

      case '?':
        return 1;
        break;

      default:
        std::cout << "Invalid argument" << std::endl;
        return 1;
    }
  }

  _device_map = get_device_map();

  if(_curve_name.empty()) {
    std::cout << "--curve required" << std::endl;
    return 1;
  }

  if(_data_dir.empty()) {
    std::cout << "--data required" << std::endl;
    return 1;
  }

  if(_curve_name != "ecp131" && _curve_name != "ecp79") {
    std::cout << "Invalid curve " << _curve_name << std::endl;
    return 1;
  }

  if(_curve_bits.find(_curve_name) == _curve_bits.end()) {
    std::cout << "Invalid curve '" << _curve_name << "'" << std::endl;
    return 1;
  }

  if(_curve_name == "ecp131") {
    _dpbits = 32;
  } else if(_curve_name == "ecp79") {
    _dpbits = 24;
  }

  if(_use_mpi && gpu_flag) {
    std::cout << "-g flag incompable when using MPI" << std::endl;
    return 1;
  }

  ecc::set_curve(_curve_name);

  // Check device ID
  int device_count = 0;
  HIP_CALL(hipGetDeviceCount(&device_count));

  if(device_count == 0) {
    std::cout << "No GPUs available" << std::endl;
    return 1;
  }


#ifdef BUILD_MPI

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

  // Use device UUID for data file
  _data_file = fmt::format("GPU-{}.dat", _device_map[_hip_device].uuid);

  // Set interrupt handler
  set_signal_handler(signal_handler);

#ifdef BUILD_HTTP
  // Start point reporter thread
  std::thread upload_thread;

  // If using MPI, only rank 0 reports results
  if(_use_upload && (_use_mpi == false || (_use_mpi == true && _world_rank == 0))) {
    upload_thread = std::thread(upload_thread_function);
  }
#endif

#ifdef BUILD_MPI
  // Thread for receiving MPI messages
  std::thread mpi_thread;
  if(_use_mpi == true && _world_rank == 0) {
    mpi_thread = std::thread(mpi_recv_thread_function);
  }
#endif

  // Run main loop
  main_loop();

#ifdef BUILD_MPI
  if(_use_mpi == true) {
    // Wait for all MPI processes to finish
    // TODO: This will wait for ALL MPI processes. If any died, this will block forever
    std::cout << "Waiting for MPI processes to finish..." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    // Stop msg receive thread
    _mpi_thread_running = false;
  }
#endif

#ifdef BUILD_HTTP
  //Wait for point thread to finish
  if(_use_upload && (_use_mpi == false || (_use_mpi == true && _world_rank == 0))) {
      upload_thread.join();
  }
#endif

#ifdef BUILD_MPI
  if(_use_mpi == true && _world_rank == 0) {
    mpi_thread.join();
  }

  // Cleanup MPI
  if(_use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}