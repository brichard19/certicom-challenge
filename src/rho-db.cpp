
#include <string.h>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <assert.h>
#include "ec_rho.h"
#include <filesystem>
#include "util.h"
#include "fmt/format.h"
#include "signal_handler.h"
#include "rocksdb/db.h"
#include "rocksdb/write_batch.h"

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

struct EncodedDP{ 
  uint8_t data[ENCODED_DP_SIZE];
};


struct JobStats {
    uint64_t num_dps;
    uint64_t total_points;
    char curve[16];
};

// The collision database is like a hashtable for distinguished points.
// Internally it stores the distingusihed points data in a list, and uses a
// set to index the list i.e. the set contains the truncated X coordinate
// and list index, and the list contains the data
class CollisionDatabase {

private:
    struct DPKey {
        uint8_t data[X_TRUNC_LEN];
    };

    struct DPData {
        uint8_t data[ENCODED_DP_SIZE - X_TRUNC_LEN];
    };
        
    rocksdb::DB* _db;

    std::vector<DPData> extract_data(const std::vector<EncodedDP> dps)
    {
        std::vector<DPData> data;

        for(auto dp : dps) {
            DPData d;

            memcpy(d.data, &dp.data[X_TRUNC_LEN], sizeof(d.data));
            data.push_back(d);
        }

        return data;
    }

    std::vector<DPKey> extract_keys(const std::vector<EncodedDP> dps)
    {
        std::vector<DPKey> keys;

        for(auto dp : dps) {
            DPKey k;

            memcpy(k.data, dp.data, sizeof(k.data));
            keys.push_back(k);
        }

        return keys;
    }

    EncodedDP construct_encoded_dp(const DPKey& key, const DPData& data)
    {
        EncodedDP encoded;

        memcpy(encoded.data, key.data, sizeof(DPKey::data));
        memcpy(encoded.data + sizeof(DPKey::data), data.data, sizeof(DPData::data));

        return encoded;
    }

public:
    CollisionDatabase(const std::string db_dir)
    {
        std::filesystem::create_directories(db_dir); 
        rocksdb::Options options;
        options.create_if_missing = true;

        rocksdb::Status status = rocksdb::DB::Open(options, db_dir + "/mydb", &_db);
        if(!status.ok()) {
            throw std::string("Error opening database");
        }
    }

    ~CollisionDatabase()
    {
        delete _db;
    }

    std::vector<std::pair<DistinguishedPoint, DistinguishedPoint>> insert(const std::vector<EncodedDP> dps)
    {
        std::vector<std::pair<DistinguishedPoint, DistinguishedPoint>> dup_dps;
   
        // Separate the Distinguished Points into key/value pairs
        std::vector<DPData> dp_data = extract_data(dps);
        std::vector<DPKey> dp_keys = extract_keys(dps);

        assert(dp_data.size() == dp_keys.size());

        rocksdb::WriteBatch batch;

        // TODO: This might not be optimal, but it works.
        for(int i = 0; i < dp_data.size(); i++) {
            rocksdb::Slice key((const char*)dp_keys[i].data, sizeof(dp_keys[i].data));
            rocksdb::Slice value((const char*)dp_data[i].data, sizeof(dp_data[i].data));

            // Check if exists
            std::string existing;
            rocksdb::Status status = _db->Get(rocksdb::ReadOptions(), key, &existing);

            if(status.IsNotFound()) {
                batch.Put(key, value);
            } else if(status.ok() == false) {
                throw std::runtime_error("rocksdb error: " + status.ToString());
            } else {
                // Found possible collision
                DPData duplicate;
                memcpy(duplicate.data, existing.data(), sizeof(duplicate.data));

                EncodedDP edp1 = construct_encoded_dp(dp_keys[i], duplicate);
                EncodedDP edp2 = construct_encoded_dp(dp_keys[i], dp_data[i]);

                DistinguishedPoint dp1 = decode_dp(edp1.data, DP_BITS);
                DistinguishedPoint dp2 = decode_dp(edp2.data, DP_BITS);
                dup_dps.push_back(std::pair<DistinguishedPoint, DistinguishedPoint>(dp1, dp2));
            }
        }

        rocksdb::Status status = _db->Write(rocksdb::WriteOptions(), &batch);
        if(status.ok() == false) {
            throw std::runtime_error("rocksdb error: " + status.ToString());
        }

        return dup_dps;
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

    memcpy(&len, &encoded.data[LEN_OFFSET], 5);

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
    of << to_str(p1.p.x) << " " << to_str(p1.p.y) << " " << to_str(p1.a) << " " << p1.length << " " << to_str(p2.a) << " " << p2.length << std::endl;
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

void main_loop()
{
  std::cout << "Monitoring " << _data_dir << " for files..." << std::endl;

  CollisionDatabase coll_db(_db_dir);

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
        decode_dp(dp.data, header.dbits);

        count += extract_length(dp);
      }

      // Keep count of total points and distinguished points
      _stats.total_points += count;
      _stats.num_dps += dps.size();

      double t0 = util::get_time();
      auto collisions = coll_db.insert(dps);
      double t1 = util::get_time();
      double time_sec = (t1 - t0);

      std::string status = fmt::format("{:<15} | DPs: {:>8} | Ps: {:>12} | Time: {:>6.3f}s", file_path, dps.size(), count, time_sec);
      std::cout << status << std::endl;

      // Process any collisions.
      if(collisions.size() > 0) {
        std::cout << "=========================================================================" << std::endl;
        for(auto coll : collisions) {
            process_collision(coll.first, coll.second);
        }
        std::cout << "=========================================================================" << std::endl;
      }

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