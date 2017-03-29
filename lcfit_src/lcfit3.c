#include "lcfit3.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "lcfit3_nlopt.h"

void lcfit3_print_array(const char* name, const size_t n, const double* x)
{
    char* sep = "";
    fprintf(stderr, "%s = { ", name);
    for (size_t i = 0; i < n; ++i) {
        fprintf(stderr, "%s%g", sep, x[i]);
        sep = ", ";
    }
    fprintf(stderr, " }\n");
}

double lcfit3_var_r(const lcfit3_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double theta_b = model->theta_b;
    const double f_1 = model->d1;

    const double r = (f_1 * (pow(theta_b, 2.0) - 1)) / ((m - c) * theta_b + m + c);

    return r;
}

double lcfit3_var_q(const lcfit3_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double theta_b = model->theta_b;

    const double q = (c - m) * theta_b - c - m;

    return q;
}

double lcfit3_var_theta(const double t, const lcfit3_bsm_t* model)
{
    const double theta_b = model->theta_b;

    const double r = lcfit3_var_r(model);
    const double theta = theta_b * exp(r * t);

    return theta;
}

void lcfit3n_gradient(const double t, const lcfit3_bsm_t* model, double* grad)
{
    const double c = model->c;
    const double m = model->m;
    const double theta_b = model->theta_b;
    const double f_1 = model->d1;

    const double r = lcfit3_var_r(model);
    const double q = lcfit3_var_q(model);
    const double theta = lcfit3_var_theta(t, model);

    //
    // This is the gradient of the normalized log-likelihood function
    // f(t) - f(0).
    //

    grad[0] = c*r*t*(theta_b - 1)/(q*theta*(1 + 1.0/theta)) + m*r*t*(theta_b - 1)/(q*theta*(-1 + 1.0/theta)) + log(1 + 1.0/theta) - log(1 + 1.0/theta_b);
    grad[1] = -c*r*t*(theta_b + 1)/(q*theta*(1 + 1.0/theta)) - m*r*t*(theta_b + 1)/(q*theta*(-1 + 1.0/theta)) + log(1 - 1/theta) - log(1 - 1/theta_b);
    grad[2] = c*((2*f_1*t*theta_b/q + r*t*(c - m)/q)/theta - 1/(theta*theta_b))/(1 + 1.0/theta) + c/(pow(theta_b, 2)*(1 + 1.0/theta_b)) + m*((2*f_1*t*theta_b/q + r*t*(c - m)/q)/theta - 1/(theta*theta_b))/(-1 + 1.0/theta) + m/(pow(theta_b, 2)*(-1 + 1.0/theta_b));
}

double lcfit3_lnl(const double t, const lcfit3_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double theta_b = model->theta_b;

    const double theta = lcfit3_var_theta(t, model);

    const double lnl = c * log(1 + 1/theta) + m * log(1 - 1/theta) - (c + m) * log(2);

    return lnl;
}

double lcfit3_norm_lnl(const double t, const lcfit3_bsm_t* model)
{
    return lcfit3_lnl(t, model) - lcfit3_lnl(0.0, model);
}

void lcfit3_evaluate_fn(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                        const size_t n, const double* t, double* lnl)
{
    for (size_t i = 0; i < n; ++i) {
        lnl[i] = lnl_fn(t[i], lnl_fn_args);
    }
}

double lcfit3_compute_weights(const size_t n, const double* lnl,
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

int lcfit3n_fit(const size_t n, const double* t, const double* lnl,
                lcfit3_bsm_t* model)
{
    double* w = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; ++i) {
        w[i] = 1.0;
    }

    int status = lcfit3n_fit_weighted(n, t, lnl, w, model);

    free(w);

    return status;
}

int lcfit3n_fit_weighted(const size_t n, const double* t, const double* lnl,
                         const double* w, lcfit3_bsm_t* model)
{
    return lcfit3n_fit_weighted_nlopt(n, t, lnl, w, model);
}

void lcfit3_normalize(const double max_lnl, const size_t n, double* lnl)
{
    for (size_t i = 0; i < n; ++i) {
        lnl[i] -= max_lnl;
    }
}

int lcfit3_fit_auto(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                    lcfit3_bsm_t* model, const double min_t, const double max_t,
                    const double alpha)
{
    const size_t n_points = 3;

    double* t = malloc(n_points * sizeof(double));
    double* lnl = malloc(n_points * sizeof(double));
    double* w = malloc(n_points * sizeof(double));

    const double max_lnl = lnl_fn(model->t0, lnl_fn_args);

    //
    // first pass
    //

    // initialize sample points

    // TODO: t0, 1/2 derivative point?, max_t
    // ...

    // here's what it looked like for lcfit2
    //lcfit2_three_points(model, lcfit2_delta(model), min_t, max_t, t);
    //t[3] = max_t;

    // evaluate, normalize, compute weights, and fit

    lcfit3_evaluate_fn(lnl_fn, lnl_fn_args, n_points, t, lnl);
    lcfit3_normalize(max_lnl, n_points, lnl);
    lcfit3_compute_weights(n_points, lnl, alpha, w);

#ifdef LCFIT3_VERBOSE
    lcfit3_print_array("t", n_points, t);
    lcfit3_print_array("w", n_points, w);
#endif /* LCFIT3_VERBOSE */

    int status = lcfit3n_fit_weighted(n_points, t, lnl, w, model);

    free(t);
    free(lnl);
    free(w);

    return status;
}
