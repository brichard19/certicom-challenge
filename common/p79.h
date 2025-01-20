#ifndef _P79_H
#define _P79_H

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

#endif