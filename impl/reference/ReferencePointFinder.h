#ifndef _REFERENCE_POINT_FINDER_H
#define _REFERENCE_POINT_FINDER_H

#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>


#include "ec_rho.h"
#include "mont.h"
#include "util.h"

class ReferencePointFinder : public DistinguishedPointFinder {

private:

  const int DISTINGUISHED_BITS = 24;
  
  const int _batch_size = 32;

  std::string _filename = "";

  std::vector<uint131_t> _a;
  std::vector<ecc::ecpoint_t> _p;
  std::vector<RWPoint> _rw;

  std::function<void(const std::vector<DistinguishedPoint>&)> _callback;

  volatile bool _running = false;



  void load()
  {
    _a.clear();
    _b.clear();
    _p.clear();

    std::ifstream ifs(_filename);

    if(!ifs.good()) {
      throw std::runtime_error("File error");
    }

    std::string line;
    while(std::getline(ifs, line)) {
      std::string a_hex;
      std::string b_hex;
      std::string x_hex;
      std::string y_hex;

      std::istringstream iss(line);
      if(!(iss >> a_hex >> b_hex >> x_hex >> y_hex)) {
        throw std::runtime_error("Invalid data");
      }

      uint131_t a = make_uint131(a_hex);
      uint131_t b = make_uint131(b_hex);
      ecc::ecpoint_t p(make_uint131(x_hex), make_uint131(y_hex));

      _a.push_back(a);
      _b.push_back(b);
      _p.push_back(p);

      assert(ecc::exists(p));
    }
  }

public:

  ReferencePointFinder()
  {
    assert(_batch_size > 0);
    assert(_batch_size == 1 || (_batch_size & (_batch_size - 1)) == 0);
  }

  size_t work_per_step()
  {
    return _batch_size;
  }

  void set_callback(std::function<void(const std::vector<DistinguishedPoint>&)> callback)
  {
    _callback = callback;
  }
  
  void init()
  {
    init("");
  }

  void init(const std::string& file)
  {
    if(!file.empty() && util::file_exists(_filename)) {
      load();
    } else {

      // Generate private keys for random starting points
      for(int i = 0; i < _batch_size; i++) {
        _a.push_back(ecc::genkey());
        _b.push_back(ecc::genkey());
      }

      // Multiply and add to get the random points 
      for(int i = 0; i < _batch_size; i++) {
        _p.push_back(ecc::mul(_a[i], ecc::g()));
        assert(ecc::exists(_p[i]));
      }
    }

    // Load the hard-coded RW points
    _rw = get_rw_points();
  }

  void step()
  {
    uint32_t dp_mask = (1 << DISTINGUISHED_BITS) - 1;

    uint32_t idx_mask = _rw.size() - 1;

    double t0 = util::get_time();

    if(_batch_size == 1) {
      for(int i = 0; i < _batch_size; i++) {

        if((_p[i].x.v[0] & dp_mask) == 0) {

          // Report distinguished point
          DistinguishedPoint dp(_a[i], _p[i]);
          std::vector<DistinguishedPoint> dps;
          dps.push_back(dp);
          _callback(dps);

          // Generate new starting point
          _a[i] = ecc::genkey();

          _p[i] = ecc::add(ecc::mul(_a[i], ecc::g()), ecc::mul(_b[i], ecc::q()));
        } else {
          uint32_t idx = _p[i].x.v[0] & idx_mask;

          _p[i] = ecc::add(_p[i], _rw[idx].p);
        }
        
        assert(ecc::exists(_p[i]));
      }
    } else {

      // Start addition
      std::vector<uint131_t> chain(_batch_size);

      uint131_t inverse = mont::to(make_uint131((uint32_t)1));

      // Compute compute 1/(Px - rx) using batch inversion.
      // Multipliy the differences together then compute
      // the inverse
      for(int i = 0; i < _batch_size; i++) {

        if((_p[i].x.v[0] & dp_mask) == 0) {
          // Report distinguished point
          DistinguishedPoint dp(_a[i], _p[i]);
          std::vector<DistinguishedPoint> dps;
          dps.push_back(dp);

          _callback(dps);

          // Generate new starting point
          _a[i] = ecc::genkey();

          _p[i] = ecc::mul(_a[i], ecc::g());
        }

        // Select R point
        int idx = _p[i].x.v[0] & idx_mask;

        // X diff 
        uint131_t diff = mont::sub(_p[i].x, _rw[idx].p.x);

        // Multiply together
        inverse = mont::mul(inverse, diff);
        chain[i] = inverse;
      }

      inverse = mont::inv(inverse);

      for(int i = _batch_size - 1; i >= 0; i--) {
        uint131_t inv_diff;
        int idx = _p[i].x.v[0] & idx_mask;

        if(i >= 1) {
          inv_diff = mont::mul(inverse, chain[i - 1]);
          uint131_t diff = mont::sub(_p[i].x, _rw[idx].p.x);
          inverse = mont::mul(inverse, diff);
        } else {
          inv_diff = inverse;
        }

        // Multiply by py - ry to get lambda
        uint131_t s = mont::mul(inv_diff, mont::sub(_p[i].y, _rw[idx].p.y));
        uint131_t s2 = mont::square(s);

        uint131_t x = mont::sub(mont::sub(s2, _p[i].x), _rw[idx].p.x);
        uint131_t y = mont::sub(mont::mul(s, mont::sub(_p[i].x, x)), _p[i].y);

        _p[i].x = x;
        _p[i].y = y;

        assert(ecc::exists(ecc::ecpoint_t(x, y)));
      }
    }
  }

  void save_progress(const std::string& filename)
  {
    std::ofstream f(filename);

    for(int i = 0; i < _p.size(); i++) {
      f << to_str(_a[i]) << " " << to_str(_p[i].x) << " " << to_str(_p[i].y) << std::endl;
    }

    f.close();
  }
};

#endif