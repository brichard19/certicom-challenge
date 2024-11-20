#ifndef _ECCP131_H
#define _ECCP131_H

#include <string>

#include "uint131.h"


namespace {

// Parameters for the ECCp-131 challenge. These are in montgomery form with
// r = 2**160.

const uint131_t _eccp131_p = {{0x194c43186b3abc0b, 0x8e1d43f293469e33, 0x4}};

const uint131_t _eccp131_a = {{0xe7f7f250cee8709a, 0xacd15fe1a8ec1522, 0x00}};

const uint131_t _eccp131_b = {{0xc85087e5ab4eca9e, 0xde124657d7ba5851, 0x2}};

const uint131_t _eccp131_n = {{0x7f7ed728f6b8e6f1, 0x8e1d43f293469e31, 0x4}};

const uint131_t _eccp131_gx = {{0xb137748018df6458, 0xb188ba8cd2a386d9, 0x00}};

const uint131_t _eccp131_gy = {{0x4b5a2a8f9b90cbef, 0x7dabb39a1bb0a5fa, 0x2}};

const uint131_t _eccp131_qx = {{0xc83af1fe332475e3, 0xe38a3357a4b0bb01, 0x2}};

const uint131_t _eccp131_qy = {{0xfe97756ed241b570, 0xb167247624e73021, 0x3}};

const uint131_t _eccp131_k = {{0xe0587d72985b105d, 0xf1fd54b0309e1ab9, 0x7cfd70cf}};

const uint131_t _eccp131_one = {{0x6e7743da32b6d0c7, 0x88c614d64c1a8f0b, 0x00}};

const uint131_t _eccp131_p_minus_2 = {{0x194c43186b3abc09, 0x8e1d43f293469e33, 0x4}};

// (P + 1) / 4
const uint131_t _eccp131_sqrt = {{0xc65310c61aceaf03, 0x238750fca4d1a78c, 0x1}};

// R
const uint131_t _eccp131_r = {{0x6e7743da32b6d0c7, 0x88c614d64c1a8f0b, 0x00}};

// R^2
const uint131_t _eccp131_r2 = {{0xf95d709f92600513, 0xf3d6fa1fb65ef639, 0x3}};

const int _eccp131_bits = 131;

const std::string _eccp131_curve_name = "ecp131";

}

namespace {

  // Parametrs for the ECCp-79 challenge. These are in montgomery form with
  // r = 2**160
const uint131_t _eccp79_p = {{0x5177412aca899cf5, 0x62ce, 0x00}};

const uint131_t _eccp79_a = {{0x732c9b460e3c3d, 0x1bb7, 0x00}};

const uint131_t _eccp79_b = {{0xc88edfd7d5b44610, 0x250c, 0x00}};

const uint131_t _eccp79_n = {{0x5177407b7258dc31, 0x62ce, 0x00}};

const uint131_t _eccp79_gx = {{0x8fa818f2b62053c8, 0x3e37, 0x00}};

const uint131_t _eccp79_gy = {{0xddab9a8daa5aa60b, 0x21a9}};

const uint131_t _eccp79_qx = {{0x96642c5fb8dbd341, 0x7a7}};

const uint131_t _eccp79_qy = {{0xbf3614d658e2931c, 0x426c}};

const uint131_t _eccp79_k = {{0x6e3655426732d0a3, 0xcafea9fd045a89b6, 0xe6bb05ec}};

const uint131_t _eccp79_one = {{0x5447aa703f6abc5f, 0x1358, 0x00}};

const uint131_t _eccp79_p_minus_2 = {{0x5177412aca899cf3, 0x62ce, 0x00}};

// (p - 5) / 8
const uint131_t _eccp79_sqrt = {{0xca2ee8255951339e, 0xc59, 0x00}};

// R^2
const uint131_t _eccp79_r2 = {{0x7b0baef57de52417, 0xe79, 0x00}};

const uint131_t _eccp79_r = {{0x5447aa703f6abc5f, 0x1358, 0x00}};

const int _eccp79_bits = 79;

const std::string _eccp79_curve_name = "ecp79";   
}


extern uint131_t _p;
extern uint131_t _a;
extern uint131_t _b;
extern uint131_t _n;
extern uint131_t _gx;
extern uint131_t _gy;
extern uint131_t _qx;
extern uint131_t _qy;
extern uint131_t _k;
extern uint131_t _one;
extern uint131_t _p_minus_2;
extern uint131_t _r;
extern uint131_t _r2;
extern int _bits;
extern int _words;
extern std::string _curve_name;

#endif