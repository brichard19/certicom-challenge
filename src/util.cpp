#include <chrono>
#include <cstdarg>
#include <fstream>
#include <random>
#include <stdexcept>
#include <thread>
#include <unistd.h>

#include "util.h"

namespace util {

double get_time()
{
  uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::system_clock::now().time_since_epoch())
                    .count();

  return (double)ms / 1000.0;
}

std::string get_date_time()
{
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);

  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  return std::string(buf);
}

std::string to_hex(const void* bytes, size_t count)
{
  std::string hex = "";

  for(int i = 0; i < count; i++) {
    char tmp[3] = {0};

    sprintf(tmp, "%.2x", ((uint8_t*)bytes)[i]);
    hex += std::string(tmp);
  }

  return hex;
}

bool file_exists(const std::string& path)
{
  std::ifstream f(path.c_str());

  return f.good();
}

// Simple timer class
Timer::Timer() {}

void Timer::start()
{
  _running = true;
  _start = util::get_time();
}

void Timer::stop()
{
  _running = false;
  _end = util::get_time();
}

double Timer::elapsed()
{
  if(_running) {
    return util::get_time() - _start;
  } else {
    return _end - _start;
  }
}

std::string get_hostname()
{
  char buf[1024];

  if(gethostname(buf, sizeof(buf))) {
    throw std::runtime_error("Error getting hostname: " + std::to_string(errno));
  }

  return std::string(buf);
}
}; // namespace util