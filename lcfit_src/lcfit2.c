#include "lcfit2.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>

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
    const double theta = exp(r * (t - t_0));
    const double v = (c - m) / theta;

    assert(isfinite(z));
    assert(isfinite(r));
    assert(isfinite(theta));
    assert(isfinite(v));

    assert(c - m != 0.0);
    assert(c + m != 0.0);
    assert(c != 0.0);
    assert(z != 0.0);
    assert(theta != 0.0);
    assert(c + m + v != 0.0);
    assert(c + m - v != 0.0);

    grad[0] = ((r*(t - t_0)/(c - m) + (t - t_0)*(z/(c + m) - z/c)/((c - m)*sqrt(z)))*v + 1/theta + 1)*c/(c + m + v) - ((r*(t - t_0)/(c - m) + (t - t_0)*(z/(c + m) - z/c)/((c - m)*sqrt(z)))*v + 1/theta - 1)*m/(c + m - v) - log(2*c + 2*m) + log(c + m + v) - 1;
    grad[1] = -((r*(t - t_0)/(c - m) - (t - t_0)*(z/(c + m) - z/m)/((c - m)*sqrt(z)))*v + 1/theta - 1)*c/(c + m + v) + ((r*(t - t_0)/(c - m) - (t - t_0)*(z/(c + m) - z/m)/((c - m)*sqrt(z)))*v + 1/theta + 1)*m/(c + m - v) - log(2*c + 2*m) + log(c + m - v) - 1;

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
    const double theta = exp(r * (t - t_0));
    const double v = (c - m) / theta;

    const double lnl = c * log(c + m + v) + m * log(c + m - v) - (c + m) * log(2 * (c + m));

    return lnl;
}

void lcfit2_rescale(const double t, const double lnl,
                    lcfit2_bsm_t* model)
{
    const double scale = lnl / lcfit2_lnl(t, model);
    model->c *= scale;
    model->m *= scale;
}

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

#ifdef LCFIT2_VERBOSE
    lcfit2_print_array("w", n, w);
#endif /* LCFIT2_VERBOSE */

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
        const double err = lcfit2_lnl(t[i], &model) - lnl[i];
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

int lcfit2_fit_weighted(const size_t n, const double* t, const double* lnl,
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

int lcfit2_fit_iterative(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                         lcfit2_bsm_t* model, double tolerance)
{
    const size_t n = 3;

    double* t = malloc(n * sizeof(double));
    double* lnl = malloc(n * sizeof(double));
    double* w = malloc(n * sizeof(double));

    const double t0 = model->t0;

    double prev_c = model->c;
    double prev_m = model->m;
    double delta_c;
    double delta_m;

    const size_t MAX_ITER = 100;
    size_t iter = 0;

    do {
        const double infl_t = lcfit2_infl_t(model);
        const double delta_t = infl_t - t0;

        t[0] = t0;
        t[1] = t0 + delta_t / 2.0;
        t[2] = infl_t;

        for (size_t i = 0; i < n; ++i) {
            lnl[i] = lnl_fn(t[i], lnl_fn_args);

            //w[i] = exp(lnl[i] - lnl[0]);
            w[i] = 1.0;
        }

#ifdef LCFIT2_VERBOSE
        lcfit2_print_array("w", n, w);
#endif /* LCFIT2_VERBOSE */

        lcfit2_rescale(t0, lnl[0], model);
        lcfit2_fit_weighted(n, t, lnl, w, model);

        delta_c = model->c - prev_c;
        delta_m = model->m - prev_m;

        prev_c = model->c;
        prev_m = model->m;

        ++iter;
    } while (fabs(delta_c / model->c) >= tolerance &&
             fabs(delta_m / model->m) >= tolerance &&
             iter < MAX_ITER);

    if (iter == MAX_ITER) {
        fprintf(stderr, "WARNING: maximum number of iterations reached in lcfit2 iterative fit\n");
    }

    free(t);
    free(lnl);
    free(w);

    return true;
}

