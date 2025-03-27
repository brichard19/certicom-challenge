#include <cassert>
#include <stdexcept>

#include "ec_rho.h"

#include "rw_p79.h"
#include "rw_p131.h"


std::vector<RWPoint> get_rw_points()
{
  std::vector<RWPoint> rw_vec;

  for(int i = 0; i < 32; i++) {
    RWPoint rw;

    if(ecc::curve_name() == "ecp131") {
      rw.a = make_uint131(_p131_a_str[i]);
      rw.b = make_uint131(_p131_b_str[i]);
      rw.p = ecc::ecpoint_t(make_uint131(_p131_x_str[i]), make_uint131(_p131_y_str[i]));
    } else if(ecc::curve_name() == "ecp79") {
      rw.a = make_uint131(_p79_a_str[i]);
      rw.b = make_uint131(_p79_b_str[i]);
      rw.p = ecc::ecpoint_t(make_uint131(_p79_x_str[i]), make_uint131(_p79_y_str[i]));
    } else {
      throw std::runtime_error("No curve selected");
    }

    // Validate
    ecc::ecpoint_t p = ecc::add(ecc::mul(rw.a, ecc::g()), ecc::mul(rw.b, ecc::q()));

    assert(ecc::is_equal(p, rw.p));

    rw_vec.push_back(rw);
  }

  return rw_vec;
}