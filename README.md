## About
Implementation of Sjostrom-Daligault non-local free-energy functional (10.1103/PhysRevLett.113.155006) in PROFESS@QE2.0 (http://www.qtp.ufl.edu/ofdft/research/computation.shtml).

PROFESS_SDF.patch is a patch file for PROFESS3.0m5B (PROFESS@QE2.0). Instructions for patching PROFESS3.0 to PROFESS3.0m5B are given at http://www.qtp.ufl.edu/ofdft/research/README.ProfAtQE.2.0.1.txt.

## Installation
1. Ensure PROFESS3.0 has been patched to PROFESS@QE2.0
2. Move to the directory `PROFESS3.0m5B/Source`
3. Run `patch -p1 < /path/to/patch/PROFESS_SDF.patch`
4. Make and install PROFESS@QE2.0 as usual

## Usage
To use the Sjostrom-Daligault functional simply set `kinetic SDF` in the input file.

The individual Fortran 90 source files are also provided.
