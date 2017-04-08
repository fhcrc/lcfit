/**
 * \file lcfit3_nlopt.c
 * \brief Implementation of lcfit3 optimization using NLopt.
 */

#include "lcfit3_nlopt.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>

#include <nlopt.h>

#include "lcfit3.h"

static const size_t MAX_ITERATIONS = 1000;

void lcfit3_print_state_nlopt(double sum_sq_err, const double* x, const double* grad)
{
    size_t iter = 0;

    fprintf(stderr, "N[%4zu] rsse = %.3f", iter, sqrt(sum_sq_err));
    fprintf(stderr, ", model = { %.3f, %.3f, %.3f }",
            x[0], x[1], x[2]);
    if (grad) {
        fprintf(stderr, ", grad = { %.6f, %.6f, %.6f }",
                grad[0], grad[1], grad[2]);
    }
    fprintf(stderr, "\n");
}

/** NLopt objective function and its gradient.
 *
 * \param[in]  p     Number of model parameters.
 * \param[in]  x     Model parameters to evaluate.
 * \param[out] grad  Gradient of the objective function at \c x.
 * \param[in]  data  Observed log-likelihood data to fit.
 *
 * \return Sum of squared error from observed log-likelihoods.
 */
double lcfit3_opt_fdf_nlopt(unsigned p, const double* x, double* grad, void* data)
{
    lcfit3_fit_data* d = (lcfit3_fit_data*) data;

    const size_t n = d->n;
    const double* t = d->t;
    const double* lnl = d->lnl;
    const double* w = d->w;

    lcfit3_bsm_t model = { x[0], x[1], x[2], d->d1, d->d2 };

    double sum_sq_err = 0.0;

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
    }

    double grad_i[3];

    for (size_t i = 0; i < n; ++i) {
        const double err = lnl[i] - lcfit3_lnl(t[i], &model);

        sum_sq_err += w[i] * pow(err, 2.0);

        if (grad) {
            lcfit3_gradient(t[i], &model, grad_i);

            grad[0] -= 2 * w[i] * err * grad_i[0];
            grad[1] -= 2 * w[i] * err * grad_i[1];
            grad[2] -= 2 * w[i] * err * grad_i[2];
        }
    }

#ifdef LCFIT3_VERBOSE
    lcfit3_print_state_nlopt(sum_sq_err, x, grad);
#endif /* LCFIT3_VERBOSE */

    return sum_sq_err;
}

/** NLopt constraint function and its gradient for enforcing that \f$c > m\f$.
 *
 * NLopt expects constraint functions of the form \f$f_c(x) \leq 0\f$,
 * so we use \f$f_c(x) = m - c\f$. That \f$c\f$ must be strictly
 * greater than \f$m\f$ is handled by the SLSQP algorithm itself, as
 * the lcfit3 log-likelihood function will return \c NaN in the case
 * where \f$c = m\f$.
 *
 * \param[in]  p     Number of model parameters.
 * \param[in]  x     Model parameters to evaluate.
 * \param[out] grad  Gradient of the constraint function at \c x.
 * \param[in]  data  Observed log-likelihood data (unused).
 *
 * \return Value of the constraint function at \c x.
 */
double lcfit3_cons_cm_nlopt(unsigned p, const double* x, double* grad, void* data)
{
    const double c = x[0];
    const double m = x[1];

    if (grad) {
        grad[0] = -1.0;
        grad[1] = 1.0;
        grad[2] = 0.0;
    }

    return m - c;
}

double lcfit3_cons_regime_3_nlopt(unsigned p, const double* x, double* grad, void* data)
{
    const double c = x[0];
    const double m = x[1];
    const double theta_b = x[2];

    if (grad) {
        grad[0] = -pow(sqrt(c) + sqrt(m), 2)/pow(c - m, 2) + (sqrt(c) + sqrt(m))/(sqrt(c)*(c - m));
        grad[1] = pow(sqrt(c) + sqrt(m), 2)/pow(c - m, 2) + (sqrt(c) + sqrt(m))/(sqrt(m)*(c - m));
        grad[2] = -1.0;
    }

    return -theta_b + pow(sqrt(c) + sqrt(m), 2)/(c - m);
}

int lcfit3_fit_weighted_nlopt(const size_t n, const double* t, const double* lnl,
                              const double* w, lcfit3_bsm_t* model)
{
    lcfit3_fit_data data = { n, t, lnl, w, model->d1, model->d2 };

    const double lower_bounds[3] = { 1.0, 1.0, 1.0 };
    const double upper_bounds[3] = { INFINITY, INFINITY, INFINITY };

    nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, 3);
    nlopt_set_min_objective(opt, lcfit3_opt_fdf_nlopt, &data);
    nlopt_set_lower_bounds(opt, lower_bounds);
    nlopt_set_upper_bounds(opt, upper_bounds);

    nlopt_add_inequality_constraint(opt, lcfit3_cons_cm_nlopt, &data, 0.0);
    nlopt_add_inequality_constraint(opt, lcfit3_cons_regime_3_nlopt, &data, 0.0);

    nlopt_set_xtol_rel(opt, sqrt(DBL_EPSILON));
    nlopt_set_maxeval(opt, MAX_ITERATIONS);

    double x[3] = { model->c, model->m, model->theta_b };
    double sum_sq_err = 0.0;

    int status = nlopt_optimize(opt, x, &sum_sq_err);

    model->c = x[0];
    model->m = x[1];
    model->theta_b = x[2];

    nlopt_destroy(opt);
    return status;
}
