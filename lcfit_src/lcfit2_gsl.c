#include "lcfit2_gsl.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "lcfit2.h"

static const size_t MAX_ITERATIONS = 1000;

int lcfit2n_opt_f(const gsl_vector* x, void* data, gsl_vector* f)
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
        //
        // We expect that the observed log-likelihoods have already
        // been normalized. The error is therefore the sum of squared
        // differences between those log-likelihoods and the
        // normalized lcfit2 log-likelihoods f(t[i]) - f(t0).
        //

        const double err = lcfit2_norm_lnl(t[i], &model) - lnl[i];
        gsl_vector_set(f, i, w[i] * err);
    }

    return GSL_SUCCESS;
}

int lcfit2n_opt_df(const gsl_vector* x, void* data, gsl_matrix* J)
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
        lcfit2n_gradient(t[i], &model, grad_i);

        gsl_matrix_set(J, i, 0, w[i] * grad_i[0]);
        gsl_matrix_set(J, i, 1, w[i] * grad_i[1]);
    }

    return GSL_SUCCESS;
}

int lcfit2n_opt_fdf(const gsl_vector* x, void* data, gsl_vector* f, gsl_matrix* J)
{
    lcfit2n_opt_f(x, data, f);
    lcfit2n_opt_df(x, data, J);

    return GSL_SUCCESS;
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

int lcfit2n_fit_weighted_gsl(const size_t n, const double* t, const double* lnl,
                             const double* w, lcfit2_bsm_t* model)
{
    double x[2] = {model->c, model->m};
    gsl_vector_const_view x_view = gsl_vector_const_view_array(x, 2);

    lcfit2_fit_data data = {n, t, lnl, w, model->t0, model->d1, model->d2};

    gsl_multifit_function_fdf fdf;

    fdf.f = &lcfit2n_opt_f;
    fdf.df = &lcfit2n_opt_df;
    fdf.fdf = &lcfit2n_opt_fdf;
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
