#include <functional>
#include <iostream>
#include <vector>

#include "ecc.h"
#include "mont.h"
#include "uint131.h"

namespace {
  bool _verbose = false;
}

bool test0()
{
  // Verify multiplying any point by order n equals point at infinity
  for(int i = 0; i < 32; i++) {
    uint131_t k = ecc::genkey();
    ecc::ecpoint_t p = ecc::mul(k, ecc::g());
    uint131_t n = ecc::n();

    ecc::ecpoint_t q = ecc::mul(n, p);

    if(ecc::is_infinity(q) == false) {
      return false;
    }
  }

  return true;
}

bool test1()
{
  // Test algebraic property (a + b) * c = ac + bc
  for(int i = 0; i < 1000; i++) {
    uint131_t a = ecc::genkey();
    uint131_t b = ecc::genkey();
    uint131_t c = ecc::genkey();


    uint131_t left = mont::mul(mont::add(a, b),c);
    uint131_t right = mont::add(mont::mul(a, c), mont::mul(b, c));

    if(_verbose) {
      std::cout << to_str(a) << std::endl;
      std::cout << to_str(b) << std::endl;
      std::cout << to_str(c) << std::endl;

      std::cout << to_str(left) << " " << to_str(right) << std::endl;
    }

    if(left != right) {
      return false;
    }
  }

  return true;
}

bool test2()
{
  // Test inverse  inv(a) * a = 1;
  for(int i = 0; i < 1000; i++) {
    uint131_t a = ecc::genkey();

    uint131_t inverse = mont::inv(a);
    
    uint131_t one = mont::mul(a, inverse);

    if(_verbose) {
      std::cout << to_str(a) << " " << to_str(inverse) << std::endl;
    }

    if(mont::mul(one, a) != a) {
      return false;
    }
  }

  return true;
}


bool test3()
{
  // Test algebraic property (a - b)(a + b) = a^2 - b^2
  for(int i = 0; i < 1000; i++) {
    uint131_t a = ecc::genkey();
    uint131_t b  = ecc::genkey();

    uint131_t ls = mont::mul(mont::sub(a, b), mont::add(a, b));
    uint131_t rs = mont::sub(mont::square(a), mont::square(b));

    if(ls != rs) {
      return false;
    }
  }

  return true;
}

bool test4()
{
  ecc::ecpoint_t g = ecc::g();
  ecc::ecpoint_t q = ecc::q();

  return ecc::exists(g) && ecc::exists(q);
}

bool test5()
{
  // Test EC point multiplication
  for(int i = 0; i < 1000; i++) {
    uint131_t k = ecc::genkey();

    ecc::ecpoint_t p = ecc::mul(k, ecc::g());

    if(!ecc::exists(p)) {
      return false;
    }
  }

  return true;
}

bool test6()
{
  ecc::ecpoint_t g = ecc::g();

  // Test point addition
  for(int i = 1; i < 100; i++) {

    uint131_t k = make_uint131(i);

    ecc::ecpoint_t p1 = ecc::mul(k, g);

    ecc::ecpoint_t p2;
    for(int j = 0; j < i; j++) {
      p2 = ecc::add(p2, g);
    }

    if(!ecc::exists(p1)) {
      std::cout << "p1 does not exist" << std::endl;
      return false;
    }
    if(!ecc::exists(p2)) {
      std::cout << "p2 does not exist" << std::endl;
      return false;
    }

    if(!(is_equal(p1, p2))) {
      return false;
    }
  }

  return true;
}

bool test7()
{
  uint131_t n = ecc::n();

  ecc::ecpoint_t p = ecc::mul(n, ecc::g());

  return ecc::is_infinity(p);
}

bool test8()
{
  ecc::ecpoint_t infinity;
  ecc::ecpoint_t g = ecc::g();

  ecc::ecpoint_t sum = ecc::add(infinity, g);

  return ecc::is_equal(sum, g);
}

bool test9()
{
  ecc::ecpoint_t infinity;

  ecc::ecpoint_t sum = ecc::dbl(infinity);

  return ecc::is_infinity(sum);
}

bool test10()
{
  for(int i = 0; i < 1000; i++) {

    // Generate quadratic residue x2
    uint131_t x = mont::square(ecc::genkey());

    // Calculate both roots
    uint131_t s1 = mont::sqrt(x);
    uint131_t s2 = mont::sub(ecc::p(), s1);

    // Square both roots
    //uint131_t y1 = mont::square(s1);
    uint131_t y1 = mont::mul(s1, s1);
    //uint131_t y2 = mont::square(s2);
    uint131_t y2 = mont::mul(s2, s2);

    // Verify they come out to x2
    if(y1 != x || y2 != x) {
      return false;
    }
  }

  return true;
}

bool test11()
{
  for(int i = 1; i < 1000; i++) {
    uint131_t x = make_uint131((uint32_t)i);

    uint131_t m = mont::to(x);
    uint131_t y = mont::from(m);

    if(x != y) {
      return false;
    }
  }

  return true;
}

int main(int argc, char**argv)
{
  std::vector<std::string> curves = {"ecp131", "ecp79"};

  for(auto curve : curves ) {
    ecc::set_curve(curve);

    std::cout << "Testing curve " << curve << std::endl;
    std::vector<std::function<bool(void)>> test_functions;

    test_functions.push_back(test0);
    test_functions.push_back(test1);
    test_functions.push_back(test2);
    test_functions.push_back(test3);
    test_functions.push_back(test4);
    test_functions.push_back(test5);
    test_functions.push_back(test6);
    test_functions.push_back(test7);
    test_functions.push_back(test8);
    test_functions.push_back(test9);
    test_functions.push_back(test10);
    test_functions.push_back(test11);

    for(int i = 0; i < test_functions.size(); i++) {
      std::cout << "Test " << (i+1) << "/" << test_functions.size() << std::endl;
    //for(auto f : test_functions) {
      if(test_functions[i]()) {
        std::cout << "    PASS" << std::endl;
      } else {
        std::cout << "    FAIL" << std::endl;
      }
    }
  }

  return 0;
}