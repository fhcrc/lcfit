#ifndef LCFIT2_NLOPT_H
#define LCFIT2_NLOPT_H

#include <stddef.h>

#include "lcfit2.h"

#ifdef __cplusplus
extern "C" {
#endif

int lcfit2_fit_weighted_nlopt(const size_t n, const double* t, const double* lnl,
                              const double* w, lcfit2_bsm_t* model);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LCFIT2_NLOPT_H */
