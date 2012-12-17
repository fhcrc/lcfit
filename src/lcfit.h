/* http://www.gnu.org/software/gsl/manual/html_node/Nonlinear-Least_002dSquares-Fitting.html */
#ifndef LCFIT_H
#define LCFIT_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif


double ll(double t, double c, double m, double r, double b);
double ml_t(const double c, const double m, const double r, const double b);
double cm_scale_factor(const double t, const double l, const double c, const double m, const double r, const double b);
int fit_ll(const size_t n, const double* t, const double* l, double* x);
int fit_ll_log(const size_t n, const double* t, const double* l, double* x, FILE* log_fp);

#ifdef __cplusplus
} // extern "C"
#endif


#endif // LCFIT_H
