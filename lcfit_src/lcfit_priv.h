/**
 * \file lcfit_priv.h
 * \brief lcfit C API - internal functions for testing
 *
 * This file provides access to internal lcfit functions for use in
 * testing.
 *
 */
#ifndef LCFIT_PRIV_H
#define LCFIT_PRIV_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

bool bracket_maximum(double (*fn)(double, void*), void* fn_args,
                     double* min_t, double* max_t);

void estimate_derivatives(double (*fn)(double, void*), void* fn_args,
                          double x, double* d1, double* d2);

double find_maximum(double (*fn)(double, void*), void* fn_args,
                    double guess, double min_t, double max_t);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // LCFIT_PRIV_H
