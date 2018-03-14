# lcfit [![build status](https://travis-ci.org/matsengrp/lcfit.svg?branch=master)](https://travis-ci.org/matsengrp/lcfit)

Likelihood curve fitting by nonlinear least squares.

## Building

### Compile-time dependencies

Building the library requires [CMake](http://www.cmake.org).

Compiling the lcfit C library requires the [GNU Scientific Library](http://www.gnu.org/software/gsl/) and [NLopt](ab-initio.mit.edu/nlopt/).
Compiling the lcfit C++ extension library, the `lcfit-compare` tool, and the test suite requires a C++11-compatible compiler.
Additionally, the `lcfit-compare` tool requires the libraries `bpp-core`, `bpp-seq`, and `bpp-phyl` from [Bio++ 2.2.0](http://biopp.univ-montp2.fr/wiki/index.php/Installation).

On Debian/Ubuntu:

```shell
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

To build the `lcfit-compare` tool required for running the example and simulations, run `make lcfit-compare`.


### Running unit tests

To build and run the test suite, run `make test`.


## Running simulations

[nestly](https://github.com/fhcrc/nestly) is used to build an extensive hierarchy of directories and configuration files to measure the behavior of `lcfit` when applied to a wide variety of data.

Running the simulations requires Python 2.7, R 3.2.4, and BppSuite 2.2.0 in addition to the dependencies required for the `lcfit-compare` tool described above.
Instructions for installing the Python and R package dependencies and running the simulations can be found in `sims/README.md`.


## Basic usage

### Fitting a model

The most straightforward way of fitting a model to an empirical log-likelihood function is via the C API function `lcfit_fit_auto()` declared in `lcfit_select.h`.
This function accepts a callback function for computing the log-likelihood at a given branch length, an initial model, and the minimum and maximum branch lengths to consider.

As an example, assume you have a function for computing log-likelihood with the following signature:

``` c++
double my_lnl(tree_t* tree, int node_id, double branch_length);
```

To adapt this function for `lcfit_fit_auto()`, you might do something like the following:

``` c++
typedef struct {
    tree_t* tree;
    int node_id;
} my_lnl_data_t;

double my_lnl_callback(double branch_length, void* data)
{
    my_lnl_data_t* lnl_data = (my_lnl_data_t*) data;
    return my_lnl(lnl_data->tree, lnl_data->node_id, branch_length);
}
```

You could then use the callback and data struct like this:

``` c++
my_lnl_data_t my_lnl_data;

// populate my_lnl_data struct
my_lnl_data.tree = ...
my_lnl_data.node_id = ...

// choose initial model
bsm_t lcfit_model = DEFAULT_INIT;  // DEFAULT_INIT defined in lcfit.c

// set bounds
double min_t = ...
double max_t = ...

// fit model and estimate maximum-likelihood branch length
double ml_t = lcfit_fit_auto(&my_lnl_callback, &my_lnl_data, &lcfit_model, min_t, max_t);

```

In C++, callback data isn't limited to plain datatypes or C-style structs.
See the function `log_likelihood_callback()` in [`lcfit_cpp_src/lcfit_compare.cc`](lcfit_cpp_src/lcfit_compare.cc) for an example of using a C++ class in the callback function for computing log-likelihoods.


### Sampling from the estimated posterior distribution

The lcfit C++ API includes a simple rejection sampler for sampling from the unnormalized posterior given an lcfit model and an exponential prior.
It's used as follows:

``` c++
gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

bsm_t lcfit_model = ...  // fitted model
double lambda = ...  // exponential distribution rate parameter

lcfit::rejection_sampler sampler(rng, lcfit_model, lambda);

// single sample
double x = sampler.sample();

// multiple samples
std::vector<double> xs = sampler.sample_n(1000);
```

The rejection sampler class also has a few functions exposed for computing the likelihood, density, and cumulative density at a given branch length.
The density and cumulative density functions rely on GSL for numerical integration.
