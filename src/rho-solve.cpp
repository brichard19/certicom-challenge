#include "ec_rho.h"
#include <fstream>
#include <sstream>

int main(int argc, char** argv)
{

  if(argc != 2) {
    std::cout << "Usage: rho-solve <FILE>" << std::endl;
    return 1;
  }

  std::string path = std::string(argv[1]);

  std::ifstream f(path);

  if(!f.is_open()) {
    std::cerr << "Error opening " << path << " for reading" << std::endl;
    return 1;
  }

  // Curve name
  std::string curve_name;
  std::string x_str;
  std::string y_str;
  std::string a1_str;
  std::string a2_str;

  try {
    // First line: curve name
    std::getline(f, curve_name);

    // Second line: space-delimited values
    std::string line;
    std::getline(f, line);

    std::istringstream iss(line);
    if(!(iss >> x_str >> y_str >> a1_str >> a2_str)) {
      throw std::runtime_error("Invalid data format on second line");
    }
  } catch(const std::ios_base::failure& e) {
    std::cerr << "File I/O error: " << e.what() << '\n';
  } catch(const std::runtime_error& e) {
    std::cerr << "Parse error: " << e.what() << '\n';
  }

  uint131_t x = make_uint131(x_str);
  uint131_t y = make_uint131(y_str);
  uint131_t a1 = make_uint131(a1_str);
  uint131_t a2 = make_uint131(a2_str);
  ecc::ecpoint_t dp(x, y);
  ecc::set_curve(curve_name);

  RhoSolver solver(dp, a1, a2);
  solver.solve();

  return 0;
}