/**
 * \file lcfit3.h
 * \brief lcfit3 C API.
 */

#ifndef LCFIT3_H
#define LCFIT3_H

#include <stddef.h>

#include "lcfit.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    /** Number of constant sites. */
    double c;
    /** Number of mutated sites. */
    double m;
    /** Free parameter dependent on mutation rate and branch length offset. */
    double theta_b;
    /** First derivative of the log-likelihood function. */
    const double d1;
    /** Second derivative of the log-likelihood function. */
    const double d2;
} lcfit3_bsm_t;

typedef struct {
    /** Number of samples. */
    const size_t n;
    /** Sampled branch lengths. */
    const double* t;
    /** Sampled log-likelihoods. */
    const double* lnl;
    /** Sample weights. */
    const double* w;
    /** First derivative of the log-likelihood function. */
    const double d1;
    /** Second derivative of the log-likelihood function. */
    const double d2;
} lcfit3_fit_data;

/** Converts an lcfit3 model to an lcfit4 model. */
void lcfit3_to_lcfit4(const lcfit3_bsm_t* model3, bsm_t* model4);

/** Computes the gradient of the log-likelihood function at branch length \c t for a given model. */
void lcfit3_gradient(const double t, const lcfit3_bsm_t* model, double* grad);

/** Computes the log-likelihood at branch length \c t for a given model. */
double lcfit3_lnl(const double t, const lcfit3_bsm_t* model);

/** Fits a model to log-likelihood data, without weighting. */
int lcfit3_fit(const size_t n, const double* t, const double* lnl,
               lcfit3_bsm_t* model);

/** Fits a model to log-likelihood data, with weighting. */
int lcfit3_fit_weighted(const size_t n, const double* t, const double* lnl,
                        const double* w, lcfit3_bsm_t* model);

/** Fits a model to a log-likelihood function. */
int lcfit3_fit_auto(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                    lcfit3_bsm_t* model, const double min_t, const double max_t,
                    const double alpha);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LCFIT3_H */
