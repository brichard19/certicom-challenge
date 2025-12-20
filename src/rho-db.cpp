
#include <string.h>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <thread>
#include <filesystem>

#include "ec_rho.h"
#include "util.h"
#include "fmt/format.h"
#include "signal_handler.h"

#define IFSTREAM_CALL(condition)\
{\
  auto& ref = condition;\
  if(!ref) {\
    std::stringstream ss;\
    ss << "FILE error " << " " << __FILE__ << ":" << __LINE__ << ": ";\
    if(ref.bad()) {\
        ss << "I/O error";\
    } else if(ref.eof()) {\
        ss << "End of file";\
    } else {\
        ss << "Read failed";\
    }\
    ss << std::endl;\
    throw std::runtime_error(ss.str());\
  }\
}\


struct JobStats {
    uint64_t num_dps;
    uint64_t total_points;
    char curve[16];
};


class RhoDb {

private:
    static const int NUM_BUCKETS = 256;
    static const int MAX_CACHE_SIZE = 64 * 1024;

    struct DBRecord {
        uint8_t key[X_TRUNC_LEN];
        DPData data;

        bool operator<(const DBRecord& other)
        {
            return memcmp(key, other.key, sizeof(key)) < 0;
        }
    };
    
    std::function<void(DistinguishedPoint, DistinguishedPoint)> _callback;

    std::array<std::mutex, NUM_BUCKETS> _locks;
    std::array<bool, NUM_BUCKETS> _dirty;
 
    std::string _db_path;
    std::vector<DBRecord> _cache[NUM_BUCKETS];
    std::thread _coll_check_thread;
    bool _running = true;

    int get_bucket(const DBRecord& rec)
    {
        unsigned int x;
        memcpy(&x, rec.key, sizeof(x));

        return x % NUM_BUCKETS;
    }

    void flush_cache(int bucket)
    {
        std::string fname = fmt::format("{}/bucket{}.dat", _db_path, bucket);

        std::ofstream f(fname, std::ios::app);

        f.write((const char*)_cache[bucket].data(), sizeof(DBRecord) * _cache[bucket].size());
        f.close();

        _cache[bucket].clear();
    }

    void flush_cache(bool force)
    {
        for(int i = 0; i < NUM_BUCKETS; i++) {
            _locks[i].lock();
            if(force || _cache[i].size() >= MAX_CACHE_SIZE) {
                flush_cache(i);
                _dirty[i] = true;
            }
            _locks[i].unlock();
        }
    }

    EncodedDP construct_encoded_dp(const uint8_t* tx, const DPData& data)
    {
        EncodedDP encoded;

        memcpy(encoded.tx, tx, X_TRUNC_LEN);

        encoded.data = data;

        return encoded;
    }
    
    void check_for_collision(int bucket)
    {
        std::string fname = fmt::format("{}/bucket{}.dat", _db_path, bucket);
        

        std::ifstream f(fname, std::ios::binary);

        // Skip if doesn't exist
        if(!f) {
            std::cout << "Skipping " << fname << std::endl;
            return;
        }
        std::cout << "Checking " << fname << " for collisions   ";

        f.seekg(0, std::ios::end);
        size_t count = f.tellg() / sizeof(DBRecord);
        f.seekg(0);

        std::vector<DBRecord> recs(count);
        f.read((char*)recs.data(), count * sizeof(DBRecord));

        f.close();

        std::cout << count << " items" << std::endl;

        if(count >= 2) {
            // Sort the records
            std::sort(recs.begin(), recs.end());

            // Check for collision
            for(int i = 0; i < recs.size() - 1; i++) {
                if(memcmp(recs[i].key, recs[i + 1].key, sizeof(recs[i].key)) == 0) {
                    // Collision found

                    EncodedDP edp1 = construct_encoded_dp(recs[i].key, recs[i].data);
                    EncodedDP edp2 = construct_encoded_dp(recs[i+1].key, recs[i+1].data);

                    DistinguishedPoint dp1 = decode_dp(edp1, DP_BITS);
                    DistinguishedPoint dp2 = decode_dp(edp2, DP_BITS);

                    _callback(dp1, dp2);
                }
            }
        }
    }

    void thread_function()
    {
        int bucket = 0;
        double last_flush = util::get_time();

        _dirty.fill(true);

        while(_running) {
            sleep(5);

            for(int bucket = 0; _running && bucket < NUM_BUCKETS; bucket++) {
                // Flush to disk every 5 minutes
                if(util::get_time() - last_flush >= 300.0) {
                    flush_cache(true);
                    last_flush = util::get_time();
                } else {
                    flush_cache(false);
                }

                if(_dirty[bucket]) {
                    check_for_collision(bucket);
                    _dirty[bucket] = false;
                }
            }
        }

        flush_cache(true);
    }

public:
    RhoDb(const std::string db_path, std::function<void(DistinguishedPoint, DistinguishedPoint)> callback) : _db_path(db_path), _callback(callback)
    {
        std::filesystem::path p = _db_path; 
        std::filesystem::create_directories(_db_path);

        _coll_check_thread = std::thread(&RhoDb::thread_function, this);
    }

    ~RhoDb()
    {
        _running = false;
        _coll_check_thread.join();
    }

    void insert(const std::vector<EncodedDP> dps)
    {
        std::vector<DBRecord> buckets[NUM_BUCKETS];

        // Sort the points into different buckets
        for(auto& dp : dps) {
            DBRecord rec;

            memcpy(rec.key, dp.tx, sizeof(rec.key));
            rec.data = dp.data;

            buckets[get_bucket(rec)].push_back(rec);
        }

        // Write to cache
        for(int b = 0; b < NUM_BUCKETS; b++) {
            _locks[b].lock();
            _cache[b].insert(_cache[b].end(), buckets[b].begin(), buckets[b].end());
            _locks[b].unlock();
        }
    }
};

namespace {

    bool _running = true;
    std::string _data_dir = "";
    std::string _db_dir = "";
    JobStats _stats;
};

void load_stats()
{
    if(std::filesystem::exists(_db_dir + "/" + "stats.bin")) {
        std::ifstream f(_db_dir + "/" + "stats.bin", std::ios::binary);

        f.read((char*)&_stats, sizeof(_stats));
        f.close();

        std::string curve_name(_stats.curve);
        ecc::set_curve(curve_name);
    }
}

void save_stats()
{
    std::ofstream f(_db_dir + "/" + "stats.bin", std::ios::binary);

    // Save curve name
    std::string curve_name = ecc::curve_name();
    sprintf(_stats.curve, "%s", curve_name.c_str());

    f.write((char*)&_stats, sizeof(_stats));
    f.close();
}

uint64_t extract_length(EncodedDP& encoded)
{
    uint64_t len = 0;

    memcpy(&len, encoded.len, sizeof(encoded.len));

    return len;
}

void process_collision(const DistinguishedPoint p1, const DistinguishedPoint p2)
{
    std::cout << "──────── Collision ────────" << std::endl;
    std::cout << "x : " << to_str(p1.p.x) << std::endl;
    std::cout << "y : " << to_str(p1.p.y) << std::endl;
    std::cout << "a : " << to_str(p1.a) << std::endl;
    std::cout << "b : " << to_str(p2.a) << std::endl;
    std::cout << "───────────────────────────" << std::endl;

    // Write the collision to a file. To be solved with a solving tool.
    std::string file_path = _db_dir + "/" + "coll" + std::to_string(rand()) + ".txt";
    std::ofstream of(file_path);

    of << ecc::curve_name() << std::endl;
    of << to_str(p1.p.x) << " " << to_str(p1.p.y) << " " << to_str(p1.a) << " " << to_str(p2.a) << std::endl;
}

double calc_probability(uint64_t n)
{
    double exponent = (double)n * n / (2 * pow(2, (double)ecc::curve_strength()));

    return 1.0 - exp(-exponent);
}

void display_status()
{
    double prob = calc_probability(_stats.total_points);
    std::string status = fmt::format("DPs: {:>9}  |  Ps: {:>15}  |  Collision: {:>7.4f}% ({:.2e})", _stats.num_dps, _stats.total_points, prob * 100.0, prob);
    std::cout << status << std::endl;
}

void collision_callback(DistinguishedPoint p1, DistinguishedPoint p2)
{
    std::cout << "=========================================================================" << std::endl;
    process_collision(p1, p2);
    std::cout << "=========================================================================" << std::endl;
    
}

void main_loop()
{
  std::cout << "Monitoring " << _data_dir << " for files..." << std::endl;

  RhoDb coll_db(_db_dir, collision_callback);
  memset(&_stats, 0, sizeof(_stats));
  load_stats();

  std::filesystem::create_directories(_data_dir);

  display_status();

  while(_running) {

    // Get list of files
    std::vector<std::string> files;

    for(const auto& entry : std::filesystem::directory_iterator(_data_dir)) {
      if(entry.is_regular_file() && entry.path().extension() == ".dat") {
        files.push_back(entry.path().filename());
      }
    }

    if(files.size() == 0) {
      // Sleep for 3 seconds so we don't sit in a busy loop when no files are available
      sleep(3);
      continue;
    }

    std::cout << "=========================================================================" << std::endl;
    
    // Process each file
    for(const auto& file_path : files) {

      std::ifstream f(_data_dir + "/" + file_path, std::ios::in | std::ios::binary);

      if(!f.is_open()) {
        throw std::runtime_error("Unable to open file for reading");
      }

      // Read the header
      DPHeader header;
      IFSTREAM_CALL(f.read((char*)&header, sizeof(header)));
      try {
        std::string name = ecc::get_curve_by_strength(header.curve);
        ecc::set_curve(name);
      } catch(...) {
        _running = false;
        break;
      }

      // Read the data
      std::vector<EncodedDP> dps(header.count);

      IFSTREAM_CALL(f.read((char *)dps.data(), header.count * sizeof(EncodedDP)));

      // Count the length of all the walks
      uint64_t count = 0;
      for(auto dp : dps) {
      
        // Validation
        decode_dp(dp, header.dbits, true);

        count += extract_length(dp);
      }

      // Keep count of total points and distinguished points
      _stats.total_points += count;
      _stats.num_dps += dps.size();

      double t0 = util::get_time();

      coll_db.insert(dps);
      double t1 = util::get_time();
      double time_sec = (t1 - t0);

      std::string status = fmt::format("{:<15} | DPs: {:>8} | Ps: {:>12} | Time: {:>6.3f}s", file_path, dps.size(), count, time_sec);
      std::cout << status << std::endl;

      // Remove file
      std::filesystem::remove(_data_dir + "/" + file_path);
    }
    std::cout << "=========================================================================" << std::endl;
    save_stats();
    display_status();
  }
}

void signal_handler(int signal)
{
  std::cout << "Exiting..." << std::endl;
  _running = false;
}

int main(int argc, char** argv)
{

    while(true) {
        static struct option long_options[] = {
            {"input-dir", required_argument, 0, 'i'},
            {"db-dir", required_argument, 0, 'd'},
            {NULL, 0, NULL, 0}
        };

        int opt_idx = 0;
        int c = getopt_long(argc, argv, "", long_options, &opt_idx);

        if(c == -1) {
            break;
        }

        switch(c) {
            case 'i':
            _data_dir = std::string(optarg);
            break;

            case 'd':
            _db_dir = std::string(optarg);
            break;

            case '?':
            return 1;

            default:
            return 1;
        }
    }

    if(_data_dir.empty() || _db_dir.empty()) {
        std::cout << "--input-dir and --db-dir required" << std::endl;
        return 1;
    }

    set_signal_handler(signal_handler);

    main_loop();

    return 0;
}