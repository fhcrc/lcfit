# lcfit [![build status](https://travis-ci.org/matsengrp/lcfit.svg?branch=master)](https://travis-ci.org/matsengrp/lcfit)

Likelihood curve fitting by nonlinear least squares.

## Building

### Compile-time dependencies

Building the library requires [CMake](http://www.cmake.org).

Compiling the lcfit C library requires the [GNU Scientific Library](http://www.gnu.org/software/gsl/) and [NLopt](ab-initio.mit.edu/nlopt/).
Compiling the lcfit C++ extension library, the `lcfit-compare` tool, and the test suite requires a C++11-compatible compiler.
Additionally, the `lcfit-compare` tool requires the libraries `bpp-core`, `bpp-seq`, and `bpp-phyl` from [Bio++ 2.2.0](http://biopp.univ-montp2.fr/wiki/index.php/Installation).


On Debian/Ubuntu:

```
sudo apt-get install \
    libgsl0-dev \
    libnlopt-dev \
    libbpp-core-dev \
    libbpp-seq-dev \
    libbpp-phyl-dev
```


### Compiling

Run `make` to obtain static and dynamic libraries.
Run `make doc` to build documentation (requires [Doxygen](http://doxygen.org)).


### Running unit tests

To build and run the test suite, run `make test`.


## Running simulations

[nestly](https://github.com/fhcrc/nestly) is used to build an extensive hierarchy of directories and configuration files to measure the behavior of `lcfit` when applied to a wide variety of data.

Running the simulations requires Python 2.7, R 3.2.4, and BppSuite 2.2.0 in addition to the dependencies required for the `lcfit-compare` tool described above.
Instructions for installing the Python and R package dependencies and running the simulations can be found in `sims/README.md`.
