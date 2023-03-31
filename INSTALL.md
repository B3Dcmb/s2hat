# Installation

This repo contains the project S2HAT

They are built using cmake (>= 3.18). Try `cmake --help` in case of doubt.

## S2HAT

In a generic use, S2HAT requires few libraries:

- MPI_C
- MPI_Fortran
- FFTW3
- HEALPIX which must be compiled in **FORTRAN**

Ensure these libraries are available on your system. You may as weel explicit them in a file `.bash_profile_s2hat` to load before using cmake, where you explicitely give `HEALPIX_DIR`, the path to find the healpix libraries and headers (especially for the Fortran one).

To build the library and install it at a given location (prefix), execute the commands:

```
cmake -S . -B build --install-prefix <path>
cmake --build build
cmake --install build
```
