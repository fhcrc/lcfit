#ifndef LCFIT2_GSL_H
#define LCFIT2_GSL_H

#include <stddef.h>

#include "lcfit2.h"

#ifdef __cplusplus
extern "C" {
#endif

int lcfit2n_fit_weighted_gsl(const size_t n, const double* t, const double* lnl,
                             const double* w, lcfit2_bsm_t* model);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LCFIT2_GSL_H */
