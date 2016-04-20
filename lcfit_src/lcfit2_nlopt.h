/**
 * \file lcfit2_nlopt.h
 * \brief NLopt-based lcfit2 fitting routines.
 */

#ifndef LCFIT2_NLOPT_H
#define LCFIT2_NLOPT_H

#include <stddef.h>

#include "lcfit2.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Fits a model to normalized log-likelihood data using NLopt, with weighting. */
int lcfit2_fit_weighted_nlopt(const size_t n, const double* t, const double* lnl,
                              const double* w, lcfit2_bsm_t* model);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LCFIT2_NLOPT_H */
