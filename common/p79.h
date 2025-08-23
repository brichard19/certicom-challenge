#ifndef _P79_H
#define _P79_H

#include "ec_params.h"

namespace {


CurveParameters _eccp79 = {
  .p = {{0x5177412aca899cf5, 0x62ce, 0x00}},

  .a = {{0x732c9b460e3c3d, 0x1bb7, 0x00}},

  .b = {{0xc88edfd7d5b44610, 0x250c, 0x00}},

  .n = {{0x5177407b7258dc31, 0x62ce, 0x00}},

  .gx = {{0x8fa818f2b62053c8, 0x3e37, 0x00}},

  .gy = {{0xddab9a8daa5aa60b, 0x21a9}},

  .qx = {{0x96642c5fb8dbd341, 0x7a7}},

  .qy = {{0xbf3614d658e2931c, 0x426c}},

  .k = {{0x6e3655426732d0a3, 0xcafea9fd045a89b6, 0xe6bb05ec}},

  .one = {{0x5447aa703f6abc5f, 0x1358, 0x00}},

  .p_minus_2 = {{0x5177412aca899cf3, 0x62ce, 0x00}},

  // (p - 5) / 8
  .sqrt = {{0xca2ee8255951339e, 0xc59, 0x00}},
  
  .r = {{0x5447aa703f6abc5f, 0x1358, 0x00}},

  // R^2
  .r2 = {{0x7b0baef57de52417, 0xe79, 0x00}},

  .bits = 79,

  .words = (79 + 63) / 64,

  .name = "ecp79"
};

}

#endif