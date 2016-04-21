/**
 * \file lcfit2.c
 * \brief Implementation of lcfit2.
 */

#include "lcfit2.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lcfit2_gsl.h"
#include "lcfit2_nlopt.h"

void lcfit2_print_array(const char* name, const size_t n, const double* x)
{
    char* sep = "";
    fprintf(stderr, "%s = { ", name);
    for (size_t i = 0; i < n; ++i) {
        fprintf(stderr, "%s%g", sep, x[i]);
        sep = ", ";
    }
    fprintf(stderr, " }\n");
}

double lcfit2_var_z(const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double f_2 = model->d2;

    const double z = (-f_2 * c * m) / (c + m);

    return z;
}

double lcfit2_var_r(const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;

    const double z = lcfit2_var_z(model);
    const double r = 2 * sqrt(z) / (c - m);

    return r;
}

double lcfit2_var_theta_tilde(const double t, const lcfit2_bsm_t* model)
{
    const double t_0 = model->t0;

    const double r = lcfit2_var_r(model);
    const double theta_tilde = exp(r * (t - t_0));

    return theta_tilde;
}

double lcfit2_var_theta(const double t, const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double t_0 = model->t0;

    const double r = lcfit2_var_r(model);
    const double theta = ((c + m) / (c - m)) * exp(r * (t - t_0));

    return theta;
}

double lcfit2_d1f_t(const double t, const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;

    const double r = lcfit2_var_r(model);
    const double theta = lcfit2_var_theta(t, model);

    const double d1f_t = (-c * r) / (theta + 1) + (m * r) / (theta - 1);

    return d1f_t;
}

double lcfit2_d2f_t(const double t, const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;

    const double r = lcfit2_var_r(model);
    const double theta = lcfit2_var_theta(t, model);

    const double d2f_t = (c * r * r * theta) / pow(theta + 1, 2.0) - (m * r * r * theta) / pow(theta - 1, 2.0);

    return d2f_t;
}

double lcfit2_infl_t(const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double t_0 = model->t0;

    const double r = lcfit2_var_r(model);
    const double b = -t_0 + (1.0 / r) * log((c + m) / (c - m));

    const double infl_t = -b + (1.0 / r) * log((c + m + 2 * sqrt(c * m)) / (c - m));

#ifdef LCFIT2_VERBOSE
    fprintf(stderr, "d2f(infl_t) = %g\n", lcfit2_d2f_t(infl_t, model));
#endif

    assert(fabs(lcfit2_d2f_t(infl_t, model)) <= 1e-6);

    return infl_t;
}

void lcfit2_model_assert_at(const double t, const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;

    const double z = lcfit2_var_z(model);
    const double r = lcfit2_var_r(model);
    const double theta_tilde = lcfit2_var_theta_tilde(t, model);
    const double v = (c - m) / theta_tilde;

#if 0
    fprintf(stderr, "model_assert: t = %g, c = %g, m = %g, t_0 = %g, f_2 = %g, z = %g, r = %g, theta_tilde = %g, v = %g\n",
            t, c, m, t_0, f_2, z, r, theta_tilde, v);
#endif

    assert(isfinite(z));
    assert(isfinite(r));
    assert(isfinite(theta_tilde));
    assert(isfinite(v));

    // basic assumptions
    assert(c > 0.0);
    assert(m > 0.0);
    assert(c > m);

    // -f_2 > 0, c > 0, m > 0, therefore z > 0
    assert(z > 0.0);

    // z > 0, c - m > 0, therefore r > 0
    assert(r > 0.0);

    // exp(x) > 0
    assert(theta_tilde > 0.0);

    // c - m > 0, theta_tilde > 0, therefore v > 0
    assert(v > 0.0);

    // log(c + m - v) must be valid
    assert(c + m - v > 0.0);
}

void lcfit2_gradient(const double t, const lcfit2_bsm_t* model, double* grad)
{
    const double c = model->c;
    const double m = model->m;
    const double t_0 = model->t0;

    const double z = lcfit2_var_z(model);
    const double r = lcfit2_var_r(model);
    const double theta_tilde = lcfit2_var_theta_tilde(t, model);
    const double v = (c - m) / theta_tilde;

    //lcfit2_model_assert_at(t, model);

    // normalized log-likelihood gradient
    grad[0] = ((r*(t - t_0)/(c - m) + (t - t_0)*(z/(c + m) - z/c)/((c - m)*sqrt(z)))*v + 1/theta_tilde + 1)*c/(c + m + v) - ((r*(t - t_0)/(c - m) + (t - t_0)*(z/(c + m) - z/c)/((c - m)*sqrt(z)))*v + 1/theta_tilde - 1)*m/(c + m - v) - log(2*c) + log(c + m + v) - 1;
    grad[1] = -((r*(t - t_0)/(c - m) - (t - t_0)*(z/(c + m) - z/m)/((c - m)*sqrt(z)))*v + 1/theta_tilde - 1)*c/(c + m + v) + ((r*(t - t_0)/(c - m) - (t - t_0)*(z/(c + m) - z/m)/((c - m)*sqrt(z)))*v + 1/theta_tilde + 1)*m/(c + m - v) + log(c + m - v) - log(2*m) - 1;

    //assert(isfinite(grad[0]));
    //assert(isfinite(grad[1]));

    //fprintf(stderr, "grad = { %g, %g }\n", grad[0], grad[1]);
}

double lcfit2_lnl(const double t, const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;

    const double theta_tilde = lcfit2_var_theta_tilde(t, model);
    const double v = (c - m) / theta_tilde;

    //lcfit2_model_assert_at(t, model);

    const double lnl = c * log(c + m + v) + m * log(c + m - v) - (c + m) * log(2 * (c + m));

    //assert(isfinite(lnl));

    return lnl;
}

void lcfit2_evaluate_fn(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                        const size_t n, const double* t, double* lnl)
{
    for (size_t i = 0; i < n; ++i) {
        lnl[i] = lnl_fn(t[i], lnl_fn_args);
    }
}

double lcfit2_compute_weights(const size_t n, const double* lnl,
                              const double alpha, double* w)
{
    double max_lnl = -HUGE_VAL;

    for (size_t i = 0; i < n; ++i) {
        if (lnl[i] > max_lnl) {
            max_lnl = lnl[i];
        }
    }

    for (size_t i = 0; i < n; ++i) {
        w[i] = pow(exp(lnl[i] - max_lnl), alpha);
    }

    return max_lnl;
}

int lcfit2_fit(const size_t n, const double* t, const double* lnl,
               lcfit2_bsm_t* model)
{
    double* w = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; ++i) {
        w[i] = 1.0;
    }

    int status = lcfit2_fit_weighted(n, t, lnl, w, model);

    free(w);

    return status;
}

int lcfit2_fit_weighted(const size_t n, const double* t, const double* lnl,
                        const double* w, lcfit2_bsm_t* model)
{
    //return lcfit2_fit_weighted_gsl(n, t, lnl, w, model);
    return lcfit2_fit_weighted_nlopt(n, t, lnl, w, model);
}

double lcfit2_delta(const lcfit2_bsm_t* model) {
    const double delta = lcfit2_infl_t(model) - model->t0;

    return delta;
}

// t must point to an array of size at least 3
void lcfit2_three_points(const lcfit2_bsm_t* model, const double delta,
                         const double min_t, const double max_t, double* t)
{
    const double t0 = model->t0;

    t[0] = t0 - delta;
    if (t[0] < min_t) {
        t[0] = (min_t + t0) / 2.0;
    }

    t[1] = t0;

    t[2] = t0 + delta;
    if (t[2] > max_t) {
        t[2] = (t0 + max_t) / 2.0;
    }
}

double lcfit2_evaluate(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                       const double alpha, const size_t n_points,
                       const double* t, double* lnl, double* w)
{
    lcfit2_evaluate_fn(lnl_fn, lnl_fn_args, n_points, t, lnl);
    double max_lnl = lcfit2_compute_weights(n_points, lnl, alpha, w);

    return max_lnl;
}

int lcfit2_fit_auto(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                    lcfit2_bsm_t* model, const double min_t, const double max_t,
                    const double alpha)
{
    const size_t n_points = 4;

    double* t = malloc(n_points * sizeof(double));
    double* lnl = malloc(n_points * sizeof(double));
    double* w = malloc(n_points * sizeof(double));

    //
    // initialize four starting points
    //

    double delta = lcfit2_delta(model);
    lcfit2_three_points(model, delta, min_t, max_t, t);

    t[3] = max_t;

    //
    // first pass
    //

    double max_lnl;
    int status;

    max_lnl = lcfit2_evaluate(lnl_fn, lnl_fn_args, alpha, n_points, t, lnl, w);

#ifdef LCFIT2_VERBOSE
    fprintf(stderr, "initial delta = %g\n", delta);
    lcfit2_print_array("t", n_points, t);
    lcfit2_print_array("w", n_points, w);
#endif

    status = lcfit2_fit_weighted(n_points, t, lnl, w, model);

    //
    // update delta, recompute the three center points, reevaluate, and refit
    //

    delta = lcfit2_delta(model);
    lcfit2_three_points(model, delta, min_t, max_t, t);

    lcfit2_evaluate(lnl_fn, lnl_fn_args, alpha, n_points, t, lnl, w);

#ifdef LCFIT2_VERBOSE
    fprintf(stderr, "revised delta = %g\n", delta);
    lcfit2_print_array("t", n_points, t);
    lcfit2_print_array("w", n_points, w);
#endif

    status = lcfit2_fit_weighted(n_points, t, lnl, w, model);

    free(t);
    free(lnl);
    free(w);

    return status;
}
