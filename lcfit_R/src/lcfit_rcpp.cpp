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
NumericVector rcpp_bsm_rescale(const double bl, const double ll, List model)
{
    // bl - branch lengths 
    // ll - log-likelihood values
    bsm_t m = {model["c"],model["m"],model["r"],model["b"]};

    lcfit_bsm_rescale(bl, ll, &m);

    return(Rcpp::NumericVector::create(
	       _["c"]=m.c,
	       _["m"]=m.m,
	       _["r"]=m.r,
	       _["b"]=m.b));
}


// [[Rcpp::export("lcfit_fit_bsm")]]
NumericVector rcpp_fit_bsm(NumericVector bl, NumericVector ll, NumericVector w, List model)
{
    // bl - branch lengths 
    // ll - log-likelihood values
    int t_n = bl.size();
    int l_n = ll.size();
    bsm_t m = {model["c"],model["m"],model["r"],model["b"]};

    assert(t_n == l_n);
    lcfit_fit_bsm_weight(t_n, bl.begin(), ll.begin(), w.begin(), &m);

    return(Rcpp::NumericVector::create(_["c"]=m.c, _["m"]=m.m, _["r"]=m.r, _["b"]=m.b));
}

