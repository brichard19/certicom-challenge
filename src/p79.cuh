#ifndef _P79_CUH
#define _P79_CUH

#include "uint131.cuh"


// Parametrs for the ECCp-79 challenge. These are in montgomery form with
// r = 2**131
__constant__ uint131_t _p79_p = {{0x5177412aca899cf5, 0x62ce, 0x00}};

__constant__ uint131_t _p79_a = {{0x1274945cc02f5789, 0x5320, 0x0}};

__constant__ uint131_t _p79_b = {{0x967af723a2d4326f, 0x1e5d, 0x0}};

__constant__ uint131_t _p79_n = {{0x5177407b7258dc31, 0x62ce, 0x00}};

__constant__ uint131_t _p79_k = {{0x6e3655426732d0a3, 0xcafea9fd045a89b6, 0x4}};

__constant__ uint131_t _p79_one = {{0x150619962eb12300, 0x5b4b, 0x0}};

__constant__ uint64_t _p79_p_minus_2[] = {0x5177412aca899cf3, 0x62ce, 0x00};

// R^2
const uint131_t _eccp79_r2 = {{0x45cc4bd33e72e43c, 0x1b44, 0x0}};

#endif