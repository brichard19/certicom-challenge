#ifndef _P79_CUH
#define _P79_CUH

#include "uint131.cuh"

// Parametrs for the ECCp-79 challenge. These are in montgomery form with
// r = 2**160
__constant__ uint131_t _p79_p = {{0x5177412aca899cf5, 0x62ce, 0x00}};

__constant__ uint131_t _p79_a = {{0x732c9b460e3c3d, 0x1bb7, 0x00}};

__constant__ uint131_t _p79_b = {{0xc88edfd7d5b44610, 0x250c, 0x00}};

__constant__ uint131_t _p79_n = {{0x5177407b7258dc31, 0x62ce, 0x00}};

__constant__ uint131_t _p79_k = {{0x6e3655426732d0a3, 0xcafea9fd045a89b6, 0xe6bb05ec}};

__constant__ uint131_t _p79_one = {{0x5447aa703f6abc5f, 0x1358, 0x00}};

__constant__ uint64_t _p79_p_minus_2[] = {0x5177412aca899cf3, 0x62ce, 0x00};

// R^2
const uint131_t _eccp79_r2 = {{0x7b0baef57de52417, 0xe79}};


#endif