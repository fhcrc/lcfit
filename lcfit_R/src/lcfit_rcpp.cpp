#include <Rcpp.h>
#include <assert.h>
#include <iostream>

#include "lcfit.h"

using namespace std;
using namespace Rcpp;


// [[Rcpp::export("lcfit_bsm_scale_factor")]]
double rcpp_bsm_scale_factor(const double t, const double l, const List model)
{
    bsm_t m = {model["c"],model["m"],model["r"],model["b"]};

    double v = lcfit_bsm_scale_factor(t, l, &m);
    return(v);

}

// [[Rcpp::export("lcfit_bsm_rescale")]]
NumericVector rcpp_bsm_rescale(const double t, const double l, List model)
{
    bsm_t m = {model["c"],model["m"],model["r"],model["b"]};

    lcfit_bsm_rescale(t, l, &m);

    return(Rcpp::NumericVector::create(
	       _["c"]=m.c,
	       _["m"]=m.m,
	       _["r"]=m.r,
	       _["b"]=m.b));
}


// [[Rcpp::export("lcfit_fit_bsm")]]
NumericVector rcpp_fit_bsm(NumericVector t, NumericVector l, List model)
{
    int t_n = t.size();
    int l_n = l.size();
    bsm_t m = {model["c"],model["m"],model["r"],model["b"]};

    assert(t_n == l_n);
    lcfit_fit_bsm(t_n, t.begin(), l.begin(), &m);

    return(Rcpp::NumericVector::create(_["c"]=m.c, _["m"]=m.m, _["r"]=m.r, _["b"]=m.b));
}

