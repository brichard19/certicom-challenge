#ifndef _P131_H
#define _P131_H

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

#endif