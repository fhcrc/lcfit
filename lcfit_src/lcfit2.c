#include "lcfit2.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

const double MAX_ITERATIONS = 1000;

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

#if 0
    fprintf(stderr, "c = %g, m = %g, t_0 = %g, f_2 = %g, z = %g, r = %g, theta_tilde = %g, v = %g\n",
            c, m, t_0, f_2, z, r, theta_tilde, v);
#endif

    assert(isfinite(z));
    assert(isfinite(r));
    assert(isfinite(theta_tilde));
    assert(isfinite(v));

    assert(c - m != 0.0);
    assert(c + m != 0.0);
    assert(c != 0.0);
    assert(z != 0.0);
    assert(theta_tilde != 0.0);
    assert(c + m + v != 0.0);
    assert(c + m - v != 0.0);

    // normalized log-likelihood gradient
    grad[0] = ((r*(t - t_0)/(c - m) + (t - t_0)*(z/(c + m) - z/c)/((c - m)*sqrt(z)))*v + 1/theta_tilde + 1)*c/(c + m + v) - ((r*(t - t_0)/(c - m) + (t - t_0)*(z/(c + m) - z/c)/((c - m)*sqrt(z)))*v + 1/theta_tilde - 1)*m/(c + m - v) - log(2*c) + log(c + m + v) - 1;
    grad[1] = -((r*(t - t_0)/(c - m) - (t - t_0)*(z/(c + m) - z/m)/((c - m)*sqrt(z)))*v + 1/theta_tilde - 1)*c/(c + m + v) + ((r*(t - t_0)/(c - m) - (t - t_0)*(z/(c + m) - z/m)/((c - m)*sqrt(z)))*v + 1/theta_tilde + 1)*m/(c + m - v) + log(c + m - v) - log(2*m) - 1;

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

    const double lnl = c * log(c + m + v) + m * log(c + m - v) - (c + m) * log(2 * (c + m));

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

// TODO: make point selection more robust to ensure bounds are respected
void lcfit2_select_points(const lcfit2_bsm_t* model,
                          const double min_t, const double max_t,
                          const size_t n, double* t)
{
    assert(n >= 2);

    const double t0 = model->t0;

    const double infl_t = lcfit2_infl_t(model);
    const double delta = 0.5 * (infl_t - t0);

    t[0] = fmax(t0 - delta, min_t + 0.5 * (t0 - min_t));
    t[1] = t0;

    for (size_t i = 2; i < n; ++i) {
        t[i] = t0 + (i - 1) * delta;
    }

#ifdef LCFIT2_VERBOSE
    lcfit2_print_array("t", n, t);
#endif /* LCFIT2_VERBOSE */
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

int lcfit2_opt_f(const gsl_vector* x, void* data, gsl_vector* f)
{
    lcfit2_fit_data* d = ((lcfit2_fit_data*) data);

    const size_t n = d->n;
    const double* t = d->t;
    const double* lnl = d->lnl;
    const double* w = d->w;

    lcfit2_bsm_t model = {gsl_vector_get(x, 0),
                          gsl_vector_get(x, 1),
                          d->t0,
                          d->d1,
                          d->d2};

    for (size_t i = 0; i < n; ++i) {
        // normalized log-likelihood error
        const double err = lcfit2_lnl(t[i], &model) - lcfit2_lnl(model.t0, &model) - lnl[i];
        gsl_vector_set(f, i, w[i] * err);
    }

    return GSL_SUCCESS;
}

int lcfit2_opt_df(const gsl_vector* x, void* data, gsl_matrix* J)
{
    lcfit2_fit_data* d = ((lcfit2_fit_data*) data);

    const size_t n = d->n;
    const double* t = d->t;
    const double* w = d->w;

    lcfit2_bsm_t model = {gsl_vector_get(x, 0),
                          gsl_vector_get(x, 1),
                          d->t0,
                          d->d1,
                          d->d2};

    double grad_i[2];

    for (size_t i = 0; i < n; ++i) {
        lcfit2_gradient(t[i], &model, grad_i);

        gsl_matrix_set(J, i, 0, w[i] * grad_i[0]);
        gsl_matrix_set(J, i, 1, w[i] * grad_i[1]);
    }

    return GSL_SUCCESS;
}

int lcfit2_opt_fdf(const gsl_vector* x, void* data, gsl_vector* f, gsl_matrix* J)
{
    lcfit2_opt_f(x, data, f);
    lcfit2_opt_df(x, data, J);

    return GSL_SUCCESS;
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

void lcfit2_print_state_gsl(size_t iter, const gsl_multifit_fdfsolver* s)
{
    gsl_vector* grad = gsl_vector_alloc(2);
    gsl_multifit_gradient(s->J, s->f, grad);

    fprintf(stderr, "G[%4zu] rsse = %.3f", iter, gsl_blas_dnrm2(s->f));
    fprintf(stderr, ", model = { %.3f, %.3f }",
            gsl_vector_get(s->x, 0),
            gsl_vector_get(s->x, 1));
    fprintf(stderr, ", grad = { %.6f, %.6f }",
            gsl_vector_get(grad, 0),
            gsl_vector_get(grad, 1));
    fprintf(stderr, "\n");

    gsl_vector_free(grad);
}

int lcfit2_fit_weighted_gsl(const size_t n, const double* t, const double* lnl,
                            const double* w, lcfit2_bsm_t* model)
{
    double x[2] = {model->c, model->m};
    gsl_vector_const_view x_view = gsl_vector_const_view_array(x, 2);

    lcfit2_fit_data data = {n, t, lnl, w, model->t0, model->d1, model->d2};

    gsl_multifit_function_fdf fdf;

    fdf.f = &lcfit2_opt_f;
    fdf.df = &lcfit2_opt_df;
    fdf.fdf = &lcfit2_opt_fdf;
    fdf.n = n;
    fdf.p = 2;
    fdf.params = &data;

    gsl_multifit_fdfsolver* s =
            gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, n, 2);

    gsl_multifit_fdfsolver_set(s, &fdf, &x_view.vector);

#ifdef LCFIT2_VERBOSE
    lcfit2_print_state_gsl(0, s);
#endif /* LCFIT2_VERBOSE */

    int status = GSL_CONTINUE;
    size_t iter = 0;

    while (status == GSL_CONTINUE && iter < MAX_ITERATIONS) {
        status = gsl_multifit_fdfsolver_iterate(s);
        ++iter;

#ifdef LCFIT2_VERBOSE
        lcfit2_print_state_gsl(iter, s);
#endif /* LCFIT2_VERBOSE */

        if (status) {
            break;
        }

        gsl_vector* grad = gsl_vector_alloc(2);
        gsl_multifit_gradient(s->J, s->f, grad);

        status = gsl_multifit_test_gradient(grad, sqrt(DBL_EPSILON));

        gsl_vector_free(grad);
    }

#ifdef LCFIT2_VERBOSE
    fprintf(stderr, "[G] status = %s (%d)   iterations %zu\n",
            gsl_strerror(status), status, iter);
#endif /* LCFIT2_VERBOSE */

    model->c = gsl_vector_get(s->x, 0);
    model->m = gsl_vector_get(s->x, 1);

    gsl_multifit_fdfsolver_free(s);

    return status;
}

int lcfit2_fit_weighted(const size_t n, const double* t, const double* lnl,
                        const double* w, lcfit2_bsm_t* model)
{
    return lcfit2_fit_weighted_gsl(n, t, lnl, w, model);
}

int lcfit2_fit_iterative(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                         lcfit2_bsm_t* model, const double min_t, const double max_t,
                         const size_t n_points, const double alpha, const size_t n_passes)
{
    const double t0 = model->t0;

    double* t = malloc(n_points * sizeof(double));
    double* lnl = malloc(n_points * sizeof(double));
    double* w = malloc(n_points * sizeof(double));

    double max_lnl;
    int status;

    for (size_t i = 0; i < n_passes; ++i) {
        lcfit2_select_points(model, min_t, max_t, n_points, t);
        lcfit2_evaluate_fn(lnl_fn, lnl_fn_args, n_points, t, lnl);
        max_lnl = lcfit2_compute_weights(n_points, lnl, alpha, w);

        // rescaling is disabled for fitting the normalized log-likelihood
        //lcfit2_rescale(t0, max_lnl, model);
        status = lcfit2_fit_weighted(n_points, t, lnl, w, model);

        // TODO: handle status codes
    }

    free(t);
    free(lnl);
    free(w);

    return status;
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

    size_t n_points = 3;

    double* t = malloc(n_points * sizeof(double));
    double* lnl = malloc(n_points * sizeof(double));
    double* w = malloc(n_points * sizeof(double));

    //
    // initialize three starting points
    //

    double delta = lcfit2_delta(model);
    lcfit2_three_points(model, delta, min_t, max_t, t);

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
