#include <Rcpp.h>
#include <assert.h>
#include <iostream>
#include <gsl/gsl_rng.h>

#include "lcfit.h"
#include "lcfit_select.h"
#include "lcfit_rejection_sampler.h"

using namespace std;
using namespace Rcpp;

// rebuild the RcppExports.cpp and RcppExports.R

// The compileAttributes function scans the source files within a package for export
// attributes and generates code as required. For example, executing this from within the
// package working directory:
//
//     compileAttributes(verbose=T)
//
// Results in the generation of the following two source files:
//
// • src/RcppExports.cpp – The extern "C" wrappers required to call exported C++
// functions within the package.
// • R/RcppExports.R – The .Call wrappers required to call the extern "C" functions
// defined in RcppExports.cpp.
//
// You should re-run compileAttributes whenever functions are added, removed, or have
// their signatures changed.

// [[Rcpp::export("lcfit_bsm_scale_factor")]]
double rcpp_bsm_scale_factor(const double t, const double l, const List model)
{
    bsm_t m = {model["c"], model["m"], model["r"], model["b"]};

    double v = lcfit_bsm_scale_factor(t, l, &m);
    return(v);

}

// [[Rcpp::export("lcfit_bsm_rescale")]]
NumericVector rcpp_bsm_rescale(const double bl, const double ll, List model)
{
    // bl - branch lengths
    // ll - log-likelihood values
    bsm_t m = {model["c"], model["m"], model["r"], model["b"]};

    lcfit_bsm_rescale(bl, ll, &m);

    return(Rcpp::NumericVector::create(
        _["c"]=m.c,
        _["m"]=m.m,
        _["r"]=m.r,
        _["b"]=m.b));
}


// [[Rcpp::export("lcfit_fit_bsm")]]
NumericVector rcpp_fit_bsm(NumericVector bl, NumericVector ll, NumericVector w, List model, int max_iter)
{
    // bl - branch lengths
    // ll - log-likelihood values
    int t_n = bl.size();
    int l_n = ll.size();
    bsm_t m = {model["c"], model["m"], model["r"], model["b"]};
    int status;

    assert(t_n == l_n);
    status = lcfit_fit_bsm_weight(t_n, bl.begin(), ll.begin(), w.begin(), &m, max_iter);

    return(Rcpp::NumericVector::create(_["c"]=m.c, _["m"]=m.m, _["r"]=m.r, _["b"]=m.b, _["status"]=status));
}

double rcpp_log_like_function(double t, void* data)
{
    Function* fn = static_cast<Function*>(data);
    return as<double>((*fn)(t));
}

// [[Rcpp::export("lcfit_fit_bsm_iter")]]
NumericVector rcpp_fit_bsm_iter(Function fn, NumericVector bl, double tolerance, List model)
{
    log_like_function_t ll_fn = {rcpp_log_like_function, &fn};
    // bl - branch lengths
    int t_n = bl.size();
    bsm_t m = {model["c"], model["m"], model["r"], model["b"]};
    bool success;

    estimate_ml_t(&ll_fn, bl.begin(), t_n, tolerance, &m, &success);

    return(Rcpp::NumericVector::create(_["c"]=m.c, _["m"]=m.m, _["r"]=m.r, _["b"]=m.b, _["success"]=success));
}

// [[Rcpp::export("lcfit_bsm_sample")]]
NumericVector rcpp_bsm_sample(List model, double lambda, int n_samples)
{
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    bsm_t m = {model["c"], model["m"], model["r"], model["b"]};

    lcfit::rejection_sampler sampler(rng, m, lambda);

    Rcpp::NumericVector samples(n_samples);

    for (int i = 0; i < n_samples; ++i) {
        samples[i] = sampler.sample();
    }

    gsl_rng_free(rng);

    return samples;
}

// exposing enums outside a class definition does not work.
// this workaround doesn't seem effective either
// http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-August/006293.html
class Dummy {
    int x;
    int get_x() {return x;}
};

RCPP_MODULE(EnumMod) {
      class_<Dummy>("EnumMod")
      .default_constructor()
      ;

      enum_<lcfit_status, Dummy>("EnumType")
          .value("LCFIT_SUCCESS", LCFIT_SUCCESS)
          .value("LCFIT_MAXITER", LCFIT_MAXITER)
          .value("LCFIT_ERROR", LCFIT_ERROR)
          ;
}
