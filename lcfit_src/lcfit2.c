#include "lcfit2.h"

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

const double MAX_ITERATIONS = 1000;


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

void lcfit2_iterative_fit(double (*lnl_fn)(double, void*),
                          void* lnl_args,
                          lcfit2_bsm_t* model,
                          const double min_t,
                          const double max_t,
                          const double tolerance,
                          bool* success)
{

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

#ifdef VERBOSE
    lcfit2_print_state_gsl(0, s);
#endif /* VERBOSE */

    int status = GSL_CONTINUE;
    size_t iter = 0;

    while (status == GSL_CONTINUE && iter < MAX_ITERATIONS) {
        status = gsl_multifit_fdfsolver_iterate(s);
        ++iter;

#ifdef VERBOSE
        lcfit2_print_state_gsl(iter, s);
#endif /* VERBOSE */

        if (status) {
            break;
        }

        status = gsl_multifit_test_delta(s->dx, s->x, 0.0, 1e-4);
    }

#ifdef VERBOSE
    fprintf(stderr, "[G] status = %s (%d)   iterations %zu\n",
            gsl_strerror(status), status, iter);
#endif /* VERBOSE */

    model->c = gsl_vector_get(s->x, 0);
    model->m = gsl_vector_get(s->x, 1);

    gsl_multifit_fdfsolver_free(s);

    return status;
}
