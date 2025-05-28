# TwistedLattice

[![Build Status](https://github.com/andrew0796/TwistedLattice.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrew0796/TwistedLattice.jl/actions/workflows/CI.yml?query=branch%3Amain)

Code for generating and analyzing fractional instantons in 4 dimensional gauge theories

## Installation
run `using Pkg; Pkg.add(path="$INSTALLATION_DIR/TwistedLattice")`

## Usage
This package defines a new type `Lattice`, which represents a four dimensional lattice gauge theory configuration with gauge group $SU(N)$, where $N$ is specified at construction, along with the size of the lattice and the twists. The lattice links can be accessed by usual indexing: to get the link at site `(i1,i2,i3,i4)` in the `mu` direction from lattice `L` you do `L[i1,i2,i3,i4,mu]`. The indexing respects the lattice periodicity in the sense that `L[i1+L1,i2,i3,i4,mu]` will return the same thing as `L[i1,i2,i3,i4,mu]` when `L1` is the length of the first direction.

Lattices are stored in [hdf5](https://support.hdfgroup.org/documentation/) files and should be read using the code here, defined in [itools.jl](https://github.com/andrew0796/TwistedLattice.jl/blob/main/src/iotools.jl). If you wish to read lattices in using a different HDF5 library, eg in python, you should be careful about the ordering of things. In particular, using `h5py` lattice sites are stored in opposite order.

Minimization is done using cooling or heat bath techniques. See the [examples](https://github.com/andrew0796/TwistedLattice.jl/tree/main/examples) and documentation therein for more details to get a sense of how to run things.

## Acknowledgements
Implementations of cooling, overrelaxation, and heat bath borrowed from https://github.com/claudio-bonati/yang-mills/blob/master/lib/sun_upd.c
Code structure heavily influenced by https://github.com/emilyzinnia/ClassicalSpinMC.jl/tree/main