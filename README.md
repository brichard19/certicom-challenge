# Certicom ECC Challenge

  

My implementaion for solving the Certicom challenge [(PDF)](https://github.com/brichard19/certicom-challenge/blob/main/docs/challenge-2009.pdf  "Certicom challenge") using GPUs.

  

# Supported curves
- P131
- P79

The target is the 131-bit prime curve. The 79-bit curve is included for testing purposes.

  

# Implementation

- Parallel Pollards rho algorithm with R-adding walks

- Distinguished points are compressed and stored in binary files
- Optional upload to central server

- All math is done in [Montgomery form](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication  "Montgomery form") with R = 2^160

- Written using HIP, so can target AMD and Nvidia GPUs

- Currently only builds on Linux

- Supports MPI for multi-GPU / multi-node parallelism

# Dependencies

### Required
* ROCm HIP SDK
### Optional 
* CUDA SDK for Nvidia GPU target
* `libopenmpi` for MPI support
* `libcurl` for server upload 
  
# Building

  Build target controlled by `TARGET_PLATFORM`. Default is "amd nvidia". 

`make TARGET_PLATFORMS=amd`or `make TARGET_PLATFORMS=nvidia`

Outputs:

`benchmark-{amd, nvidia}`

`rho-{amd, nvidia}`

# Running

## Benchmark

`benchmark- --curve ecp131`

  

## Rho program

  

`rho --curve ecp131 --data-dir `

  

### Rho program with HTTP upload

  

`rho --curve ecp131 --upload-server <hostname>`

  

### Rho program with MPI

  

`mpirun -n 5 rho --mpi --curve ecp131`