/**
 * \file lcfit2.h
 * \brief lcfit2 C API.
 */

#ifndef LCFIT2_H
#define LCFIT2_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** lcfit2 model parameters. */
typedef struct {
    /** Number of constant sites. */
    double c;
    /** Number of mutated sites. */
    double m;
    /** Maximum-likelihood branch length. */
    const double t0;
    /** First derivative of the log-likelihood function. */
    const double d1;
    /** Second derivative of the log-likelihood function. */
    const double d2;
} lcfit2_bsm_t;

typedef struct {
    /** Number of samples. */
    const size_t n;
    /** Sampled branch lengths. */
    const double* t;
    /** Sampled log-likelihoods. */
    const double* lnl;
    /** Sample weights. */
    const double* w;
    /** Maximum-likelihood branch length. */
    const double t0;
    /** First derivative of the log-likelihood function. */
    const double d1;
    /** Second derivative of the log-likelihood function. */
    const double d2;
} lcfit2_fit_data;

/** Computes the inflection point of the log-likelihood function for a given model. */
double lcfit2_infl_t(const lcfit2_bsm_t* model);

/** Computes the gradient of the log-likelihood function at branch length \c t for a given model. */
void lcfit2_gradient(const double t, const lcfit2_bsm_t* model, double* grad);

/** Computes the unnormalized log-likelihood at branch length \c t for a given model. */
double lcfit2_lnl(const double t, const lcfit2_bsm_t* model);

/** Computes the normalized log-likelihood at branch length \c t for a given model. */
double lcfit2n_lnl(const double t, const lcfit2_bsm_t* model);

/** Fits a model to normalized log-likelihood data, without weighting. */
int lcfit2_fit(const size_t n, const double* t, const double* lnl,
               lcfit2_bsm_t* model);

/** Fits a model to normalized log-likelihood data, with weighting. */
int lcfit2_fit_weighted(const size_t n, const double* t, const double* lnl,
                        const double* w, lcfit2_bsm_t* model);

/** Fits a model to a normalized log-likelihood function. */
int lcfit2_fit_auto(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                    lcfit2_bsm_t* model, const double min_t, const double max_t,
                    const double alpha);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LCFIT2_H */
