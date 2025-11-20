
#include <string.h>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <getopt.h>

#include "ec_rho.h"
#include "giga/gigalist.h"
#include "giga/gigaset.h"
#include "util.h"
#include "fmt/format.h"

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
        uint64_t id;
        uint8_t data[X_TRUNC_LEN];
    };

    struct DPData {
        uint8_t data[ENCODED_DP_SIZE - X_TRUNC_LEN];
    };

    static int compare_keys(const DPKey k1, const DPKey k2)
    {
        return memcmp(k1.data, k2.data, 13);
    }

    static uint32_t hash_key(const DPKey k, int n)
    {
        return ((uint32_t*)&k.data)[0] % n;
    }

    GigaList<DPData> _list;

    GigaSet<DPKey, CollisionDatabase::compare_keys, CollisionDatabase::hash_key> _set;

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
    CollisionDatabase(const std::string dbname) : _list(dbname + "/" + "list.dat"), _set(dbname + "/" + "index.dat")
    {

    }


    std::vector<std::pair<DistinguishedPoint, DistinguishedPoint>> insert(const std::vector<EncodedDP> dps)
    {
        std::vector<std::pair<DistinguishedPoint, DistinguishedPoint>> dup_dps;
   
        // Separate the Disntingusihed Points into key/value pairs
        std::vector<DPData> dp_data = extract_data(dps);
        std::vector<DPKey> dp_keys = extract_keys(dps);

        // Append to list first
        std::vector<uint64_t> idxs = _list.append(dp_data);

        assert(idxs.size() == dp_keys.size());

        // Store the list idx with the key
        for(int i = 0; i < idxs.size(); i++) {
            dp_keys[i].id = idxs[i];
        }

        // Convert duplicates to EncodedDP and return
        auto dup_keys = _set.insert(dp_keys);

        if(dup_keys.size() > 0) {

            for(int i = 0; i < dup_keys.size(); i++) {

                DPData data1 = _list.get(dup_keys[i].first.id);
                DPData data2 = _list.get(dup_keys[i].second.id);
                 
                EncodedDP edp1 = construct_encoded_dp(dup_keys[i].first, data1);
                EncodedDP edp2 = construct_encoded_dp(dup_keys[i].second, data2);

                DistinguishedPoint dp1 = decode_dp(edp1.data, DP_BITS);
                DistinguishedPoint dp2 = decode_dp(edp2.data, DP_BITS);

                dup_dps.push_back(std::pair<DistinguishedPoint, DistinguishedPoint>(dp1, dp2));
            }

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
      if(header.curve == 131) {
        ecc::set_curve("ecp131");
      } else if(header.curve == 79) {
        ecc::set_curve("ecp79");
      } else if(header.curve == 89) {
        ecc::set_curve("ecp89");
      } else {
        std::cout << "Invalid curve" << std::endl;
        exit(1);
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

    main_loop();

    return 0;
}