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

double lcfit2_d1f_t(const double t, const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double t_0 = model->t0;
    const double f_2 = model->d2;

    const double z = (-f_2 * c * m) / (c + m);
    const double r = 2 * sqrt(z) / (c - m);
    const double theta = ((c + m) / (c - m)) * exp(r * (t - t_0));

    const double d1f_t = (-c * r) / (theta + 1) + (m * r) / (theta - 1);

    return d1f_t;
}

double lcfit2_d2f_t(const double t, const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double t_0 = model->t0;
    const double f_2 = model->d2;

    const double z = (-f_2 * c * m) / (c + m);
    const double r = 2 * sqrt(z) / (c - m);
    const double theta = ((c + m) / (c - m)) * exp(r * (t - t_0));

    const double d2f_t = (c * r * r * theta) / pow(theta + 1, 2.0) - (m * r * r * theta) / pow(theta - 1, 2.0);

    return d2f_t;
}

double lcfit2_infl_t(const lcfit2_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double t_0 = model->t0;
    const double f_2 = model->d2;

    const double z = (-f_2 * c * m) / (c + m);
    const double r = 2 * sqrt(z) / (c - m);
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
    const double t_0 = model->t0;
    const double f_2 = model->d2;

    const double z = (-f_2 * c * m) / (c + m);
    const double r = 2 * sqrt(z) / (c - m);
    const double theta_tilde = exp(r * (t - t_0));
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
    const double f_2 = model->d2;

    const double z = (-f_2 * c * m) / (c + m);
    const double r = 2 * sqrt(z) / (c - m);
    const double theta_tilde = exp(r * (t - t_0));
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
    const double t_0 = model->t0;
    const double f_2 = model->d2;

    const double z = (-f_2 * c * m) / (c + m);
    const double r = 2 * sqrt(z) / (c - m);
    const double theta_tilde = exp(r * (t - t_0));
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

void lcfit2_rescale(const double t, const double lnl,
                    lcfit2_bsm_t* model)
{
    const double scale = lnl / lcfit2_lnl(t, model);
    model->c *= scale;
    model->m *= scale;
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

static int compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da > *db) - (*da < *db);
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
        t[0] = min_t + (t0 - min_t) / 2.0;
    }

    t[1] = t0;

    t[2] = t0 + delta;
    if (t[2] > max_t) {
        t[2] = t0 + (max_t - t0) / 2.0;
    }
}

// t must point to a sorted array of size at least 2
// omega is the weight threshold
double lcfit2_select_point_left(const double delta, const double min_t, const double omega,
                                const double* t, const double* w)
{
    double new_t;

    // test using weight of leftmost point t[0]
    if (w[0] < omega) {
        // if weight < threshold, choose point halfway between t[0] and t[1]
        new_t = t[0] + (t[1] - t[0]) / 2.0;
    } else {
        // if weight >= threshold, choose point t[0] - delta or halfway
        // between t[0] and min_t if t[0] - delta is less than min_t
        new_t = t[0] - delta;
        if (new_t < min_t) {
            new_t = min_t + (t[0] - min_t) / 2.0;
        }
    }

    return new_t;
}

// t must point to a sorted array of size at least 2
// omega is the weight threshold
double lcfit2_select_point_right(const double delta, const double max_t, const double omega,
                                 const double* t, const double* w)
{

    double new_t;

    if (w[1] < omega) {
        // if weight < threshold, choose point halfway between t[0] and t[1]
        new_t = t[0] + (t[1] - t[0]) / 2.0;
    } else {
        // if weight >= threshold, choose point t[1] + delta or halfway
        // between t[1] and max_t if t[1] + delta is greater than max_t
        new_t = t[1] + delta;
        if (new_t > max_t) {
            new_t = t[1] + (max_t - t[1]) / 2.0;
        }
    }

    return new_t;
}

double lcfit2_evaluate(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                       lcfit2_bsm_t* model, const double alpha, const size_t n_points,
                       const double* t, double* lnl, double* w)
{
    lcfit2_evaluate_fn(lnl_fn, lnl_fn_args, n_points, t, lnl);
    double max_lnl = lcfit2_compute_weights(n_points, lnl, alpha, w);

    // rescaling is disabled for fitting the normalized log-likelihood
    //lcfit2_rescale(model->t0, max_lnl, model);

    return max_lnl;
}

int lcfit2_fit_iterative2(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                          lcfit2_bsm_t* model, const double min_t, const double max_t,
                          const double alpha, const double omega, const size_t n_passes)
{
    assert(n_passes >= 1);

    size_t n_points = 4;

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

    max_lnl = lcfit2_evaluate(lnl_fn, lnl_fn_args, model, alpha, n_points, t, lnl, w);

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

    lcfit2_evaluate(lnl_fn, lnl_fn_args, model, alpha, n_points, t, lnl, w);

#ifdef LCFIT2_VERBOSE
    fprintf(stderr, "revised delta = %g\n", delta);
    lcfit2_print_array("t", n_points, t);
    lcfit2_print_array("w", n_points, w);
#endif

    status = lcfit2_fit_weighted(n_points, t, lnl, w, model);


    //
    // additional passes
    //

    for (size_t i = 1; i < n_passes; ++i) {
        //
        // compute new points using the first two and last two evaluated points
        //

        double t_left = lcfit2_select_point_left(delta, min_t, omega, t, w);
        double t_right = lcfit2_select_point_right(delta, max_t, omega, t + n_points - 2, w + n_points - 2);

        t = realloc(t, (n_points + 2) * sizeof(double));
        memmove(t + 1, t, n_points * sizeof(double));
        n_points += 2;

        t[0] = t_left;
        t[n_points - 1] = t_right;

        // here we sort the points again, because t_left and t_right
        // can be to the right or left of their neighboring points. a
        // smarter thing to do would be to swap the points above if
        // the new points are out of order. Note that we don't bother
        // reordering lnl and w accordingly since they get
        // recalculated anyway in lcfit2_evaluate.

        qsort(t, n_points, sizeof(double), compare_doubles);

        //
        // reallocate lnl and w
        //

        lnl = realloc(lnl, n_points * sizeof(double));
        w = realloc(w, n_points * sizeof(double));

        //
        // refit
        //

        lcfit2_evaluate(lnl_fn, lnl_fn_args, model, alpha, n_points, t, lnl, w);

#ifdef LCFIT2_VERBOSE
        lcfit2_print_array("t", n_points, t);
        lcfit2_print_array("w", n_points, w);
#endif

        status = lcfit2_fit_weighted(n_points, t, lnl, w, model);
    }

    free(t);
    free(lnl);
    free(w);

    return status;
}
