# Certicom ECC Challenge

My implementaion for solving the Certicom challenge [(PDF)](https://github.com/brichard19/certicom-challenge/blob/main/docs/challenge-2009.pdf "Certicom challenge") using GPUs.

# Supported curves

- p131
- p79

The target is the 131-bit prime curve. The 79-bit curve is included for testing purposes.

#Implementation
- Parallel Pollards rho algorithm with R-adding walks.
- Distinguished points are uploaded to central server. Currently the server code is not implemented, so the upload gets skipped.
- All math is done in  [Montgomery form](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication "Montgomery form") with R = 2^160
- Built using HIP, so can target AMD and Nvidia GPus
- Currently builds on Linux

# Dependencies

- ROCm HIP SDK
- CUDA SDK (for Nvidia build)
- libcurl
- Optional: libopenmpi (for MPI support)


# Building

- Edit `Makefile` and set `TARGET_PLATFORMS` appropriately.

- Run `make`

Outputs:

`bin/certicom-bench-{amd, nvidia}`
`bin/certicom-rho-{amd, nvidia}`

# Running

## Benchmark
`certicom-bench-{amd, nvidia} ecp131`

## Rho program

`certicom-rho-{amd, nvidia} ecp131`
