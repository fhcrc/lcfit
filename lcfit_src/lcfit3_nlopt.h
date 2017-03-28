/**
 * \file lcfit3_nlopt.h
 * \brief NLopt-based lcfit3 fitting routines.
 */

#ifndef LCFIT3_NLOPT_H
#define LCFIT3_NLOPT_H

#include <stddef.h>

#include "lcfit3.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Fits a model to normalized log-likelihood data using NLopt, with weighting. */
int lcfit3n_fit_weighted_nlopt(const size_t n, const double* t, const double* lnl,
                               const double* w, lcfit3_bsm_t* model);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LCFIT3_NLOPT_H */
