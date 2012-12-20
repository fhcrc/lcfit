/* http://www.gnu.org/software/gsl/manual/html_node/Nonlinear-Least_002dSquares-Fitting.html */
#ifndef LCFIT_H
#define LCFIT_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif


double lcfit_bsm_log_like(double t, double c, double m, double r, double b);
double lcfit_bsm_ml_t(const double c, const double m, const double r, const double b);
double lcfit_bsm_scale_factor(const double t, const double l, const double c, const double m, const double r, const double b);
int lcfit_fit_bsm(const size_t n, const double* t, const double* l, double* x);

#ifdef __cplusplus
} // extern "C"
#endif


#endif // LCFIT_H
