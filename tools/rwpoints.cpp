#include <iostream>
#include <cassert>

#include "ec_rho.h"
#include "fmt/format.h"
#include "util.h"

void print_rw_points()
{
  std::string name = ecc::curve_name();
  std::vector<RWPoint> rw_vec;

  std::vector<uint131_t> a_vec = ecc::genkeys(32, 1234);
  std::vector<uint131_t> b_vec = ecc::genkeys(32, 5678);

  for(int i = 0; i < 32; i++) {
    RWPoint rw;

    rw.a = a_vec[i];
    rw.b = b_vec[i];

    rw.p = ecc::add(ecc::mul(rw.a, ecc::g()), ecc::mul(rw.b, ecc::q()));

    assert(ecc::exists(rw.p));

    rw_vec.push_back(rw);
  }


  std::cout << fmt::format("std::string _{}_a_str[] = {{", name) << std::endl;
  for(int i = 0; i < rw_vec.size(); i++) {
    std::cout << fmt::format("\"{}\",", to_str(rw_vec[i].a)) << std::endl;
  }
  std::cout << "};" << std::endl;

  std::cout << std::endl;
 
  std::cout << fmt::format("std::string _{}_b_str[] = {{", name) << std::endl;
  for(int i = 0; i < rw_vec.size(); i++) {
    std::cout << fmt::format("\"{}\",", to_str(rw_vec[i].b)) << std::endl;
  }
  std::cout << "};" << std::endl;
  std::cout << std::endl;
  
  std::cout << fmt::format("std::string _{}_x_str[] = {{", name) << std::endl;
  for(int i = 0; i < rw_vec.size(); i++) {
    std::cout << "\"" << to_str(rw_vec[i].p.x) << "\"," << std::endl;
  }
  std::cout << "};" << std::endl;
  std::cout << std::endl;
  
  std::cout << fmt::format("std::string _{}_y_str[] = {{", name.c_str()) << std::endl;
  for(int i = 0; i < rw_vec.size(); i++) {
    std::cout << "\"" << to_str(rw_vec[i].p.y) << "\"," << std::endl;
  }
  std::cout << "};" << std::endl;
}

void print_multiples(ecc::ecpoint_t g, std::string name)
{
  std::vector<ecc::ecpoint_t> points;

  ecc::ecpoint_t p = g;
  for(int i = 0; i < 131; i++) {
    points.push_back(p);

    p = ecc::dbl(p);
  }

  int bits = ecc::curve_strength();

  // Print x
  std::cout << fmt::format("__device__ uint131_t _{}x[{}] = {{", name.c_str(), bits) << std::endl;
  
  for(int i = 0; i < bits; i++) {
    std::cout << fmt::format("{{0x{:016x},0x{:016x},0x{:02x}}},", points[i].x.v[0], points[i].x.v[1], points[i].x.v[2]) << std::endl;
  }

  std::cout << "};" << std::endl << std::endl;

  // Print y
  std::cout << fmt::format("__device__ uint131_t _{}y[{}] = {{", name.c_str(), bits) << std::endl;
  for(int i = 0; i < bits; i++) {
    std::cout << fmt::format("{{0x{:016x},0x{:016x},0x{:02x}}},", points[i].y.v[0], points[i].y.v[1], points[i].y.v[2]) << std::endl;
  }

  std::cout << "};" << std::endl << std::endl;
}

int main(int argc, char** argv)
{
  if(argc != 2) {
    std::cout << "Usage: rwpoints [curve]" << std::endl;
    return 1;
  }
  
  std::string curve_name = std::string(argv[1]);

  ecc::set_curve(curve_name);

  std::cout << "// Paramters for curve " << ecc::curve_name() << std::endl;
  print_rw_points();

  print_multiples(ecc::g(), "g");
 
  print_multiples(ecc::q(), "q");

  return 0;
}