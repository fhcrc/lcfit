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

void lcfit2_gradient(const double t, const lcfit2_bsm_t* model, double* grad);

double lcfit2_lnl(const double t, const lcfit2_bsm_t* model);

void lcfit2_rescale(const double t, const double lnl,
                    lcfit2_bsm_t* model);

int lcfit2_fit(const size_t n, const double* t, const double* lnl,
               lcfit2_bsm_t* model);

int lcfit2_fit_weighted(const size_t n, const double* t, const double* lnl,
                        const double* w, lcfit2_bsm_t* model);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LCFIT2_H */
