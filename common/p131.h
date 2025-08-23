#ifndef _P131_H
#define _P131_H

#include "ec_params.h"

namespace {

CurveParameters _eccp131 = {
    .p = {{0x194c43186b3abc0b, 0x8e1d43f293469e33, 0x4}},
    .a = {{0xe7f7f250cee8709a, 0xacd15fe1a8ec1522, 0x00}},
    .b = {{0xc85087e5ab4eca9e, 0xde124657d7ba5851, 0x2}},
    .n = {{0x7f7ed728f6b8e6f1, 0x8e1d43f293469e31, 0x4}},

    .gx = {{0xb137748018df6458, 0xb188ba8cd2a386d9, 0x00}},
    .gy = {{0x4b5a2a8f9b90cbef, 0x7dabb39a1bb0a5fa, 0x2}},

    .qx = {{0xc83af1fe332475e3, 0xe38a3357a4b0bb01, 0x2}},

    .qy = {{0xfe97756ed241b570, 0xb167247624e73021, 0x3}},

    // k such that k * k^-1 = -1 (mod R)
    // k = (r * r_inv -1) // p
    .k = {{0xe0587d72985b105d, 0xf1fd54b0309e1ab9, 0x7cfd70cf}},

    .one = {{0x6e7743da32b6d0c7, 0x88c614d64c1a8f0b, 0x00}},

    .p_minus_2 = {{0x194c43186b3abc09, 0x8e1d43f293469e33, 0x4}},

    // (P + 1) / 4
    .sqrt = {{0xc65310c61aceaf03, 0x238750fca4d1a78c, 0x1}},

    // R mod p
    .r = {{0x6e7743da32b6d0c7, 0x88c614d64c1a8f0b, 0x00}},

    // R^2
    .r2 = {{0xf95d709f92600513, 0xf3d6fa1fb65ef639, 0x3}},

    .bits = 131,

    .words = (131 + 63) / 64,

    .name = "ecp131"
};

}
#endif