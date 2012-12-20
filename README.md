# Building

## Dependencies

Building the library requires [cmake][1]

The lcfit C-library requires the [GNU Scientific Library][2] (`libgsl0-dev` on debian).

The `lcfit-compare` tool requires a C++11-compatible compiler, and the `bpp-core`, `bpp-seq`, and `bpp-phyl` libraries from the [Bio++ suite][3] (`libbpp-{core,seq,phyl}-dev` on debian).

Running the unit tests requires a C++11-compatible compiler.

## Compiling

Run `make`

## Running unit tests

To run the test suite, run `make test`

[1]: http://www.cmake.org
[2]: http://www.gnu.org/s/gsl
[3]: http://biopp.univ-montp2.fr
