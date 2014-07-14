# `lcfit`

Likelihood curve fitting by nonlinear least squares.

# Building

## Dependencies

Building the library requires [cmake][1]

The lcfit C-library requires the [GNU Scientific Library][2] (`libgsl0-dev` on debian).

The `lcfit-compare` tool requires a C++11-compatible compiler, and the `bpp-core`, `bpp-seq`, and `bpp-phyl` libraries from the [Bio++ suite master branch ][3].

We used:

    bpp-core
    commit 07e645acf4a90a81eb555d4a6ff7fe9ae951fd68

    bpp-phyl
    commit 498fa3fa6638bdad777e420c72873d8f3d68c00a

    bpp-seq
    commit 1ba912f48ad0eb42369269dc26209d66d7e1ccd9

    bppsuite
    commit 324d6f761c28c2fa4380a7233a59696e36b6a5e0

Running the unit tests requires python.

## Compiling

Run `make` to obtain static and dynamic libraries.

Run `make doc` to build documentation.

## Running unit tests

To run the test suite, run `make test`

[1]: http://www.cmake.org
[2]: http://www.gnu.org/s/gsl
[3]: http://biopp.univ-montp2.fr

## Running simulations ##

[Nestly](https://github.com/fhcrc/nestly) is used to build an extensive hierarchy of
directories and configuration files to measure the behavior of lcfit
when applied to a wide variety of data.  The simulations can take several hours to
complete, but they can be run in parallel by passing the `-j` option
to `scons`.   This will cause `scons` to launch multiple `slurm` jobs at
once instead of running them sequentially.

    $ cd sims
	$ scons -j 10
	
