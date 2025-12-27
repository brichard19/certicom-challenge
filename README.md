# Certicom ECC Challenge

My GPU accelerated implementaion for solving the Certicom challenge [(PDF)](/docs/challenge-2009.pdf  "Certicom challenge") using GPUs.

# Supported curves
- P131
- P89
- P79

The target is the 131-bit prime curve. Smaller curves are included for testing purposes.

# Implementation
- Parallel Pollard's rho algorithm with an R-adding walk
- Distinguished points are compressed and stored in binary files
- All math is done in [Montgomery form](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication  "Montgomery form") with R = 2^160
- Uses AMD HIP, so can target both AMD and Nvidia GPUs
- Currently only builds on Linux
- Supports MPI for multi-GPU / multi-node parallelism

 
### The project consists of 3 tools

#### rho
`rho-amd`/`rho-nvidia`

  This program runs the R-added walks on the GPUs and finds distinguished points and saves them to a file.

usage:
`rho-amd --dp-bits 20 --curve ecp131 --data-dir /path/to/data/dir`

Per-GPU data is stored in `/path/to/data/dir/hostname/GPU-<serial number>`

The results are stored in `/path/to/data/dir/results`


#### rho-db
This program reads data files output from rho program and indexes the distinguished points. It continually monitors an input directory for new data.

When two identical distinguished points are found, they are written to a file.

The database strategy is similar to the one described in  [Breaking ECC2K-130](/docs/breaking-ecc2k-130.pdf  "Breaking ECC2K-130").

Usage:
`rho-db --input-dir /path/to/data/dir/results  --db-dir /path/to/database/dir`

#### solver.py

The solver program takes the two distinguished points from rho-db and solves for the private key

Usage:
`solver.py collision_file.txt`

# Dependencies

### Required

* ROCm HIP SDK

### Optional

* CUDA SDK for Nvidia GPU target
*  `libopenmpi` for MPI support

# Building

Build target controlled by `TARGET_PLATFORM`. Default is "amd nvidia".

`make TARGET_PLATFORMS=amd`or `make TARGET_PLATFORMS=nvidia`

  
Outputs:
`benchmark-{amd, nvidia}`

`rho-{amd, nvidia}`

`rho-database`