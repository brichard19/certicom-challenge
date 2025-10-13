#ifndef _UTIL_H
#define _UTIL_H

#include <stdint.h>
#include <string>
#include <vector>

namespace util {

double get_time();
std::string get_date_time();
std::string get_hostname();
std::string to_hex(const void* bytes, size_t count);

bool file_exists(const std::string& path);

class Timer {

private:
    double _start = 0.0;

    double _end = 0.0;

    bool _running = false;

public:
    Timer();

    void start();
    void stop();
    double elapsed();
};

};

#endif