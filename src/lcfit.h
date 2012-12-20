/* http://www.gnu.org/software/gsl/manual/html_node/Nonlinear-Least_002dSquares-Fitting.html */
#ifndef LCFIT_H
#define LCFIT_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Fit of the binary-symmetric model */
typedef struct {
    double c;
    double m;
    double r;
    double b;
} bsm_t;

double lcfit_bsm_log_like(double, const bsm_t*);
double lcfit_bsm_ml_t(const bsm_t*);
double lcfit_bsm_scale_factor(const double, const double, const bsm_t*);
int lcfit_fit_bsm(const size_t, const double*, const double*, bsm_t*);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // LCFIT_H
