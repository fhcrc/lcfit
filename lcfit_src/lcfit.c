/* http://www.gnu.org/software/gsl/manual/html_node/Nonlinear-Least_002dSquares-Fitting.html */
/**
 * \file lcfit.c
 * \brief lcfit C-API implementation.
 */

#include "lcfit.h"

#include <assert.h>
#include <math.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_roots.h>

#include <nlopt.h>

const static double LAMBDA = 50;

/** Minimum bound on mutation rate. */
const static double BSM_R_MIN = 1e-9;

/** Maximum bound on mutation rate. */
const static double BSM_R_MAX = 100.0;

/** Minimum bound on branch length offset. */
const static double BSM_B_MIN = 1e-12;

const bsm_t DEFAULT_INIT = {1100.0, 800.0, 2.0, 0.5};

/* calculate the sum of logs */
static double log_sum(const double x, const double y)
{
    const static double max_value = DBL_MAX;
    const double log_limit = -max_value / 100;
    const static double NATS = 400;

    const double temp = y - x;
    if(temp > NATS || x < log_limit)
        return y;
    if(temp < -NATS || y < log_limit)
        return x;
    if(temp < 0)
        return x + log1p(exp(temp));
    return y + log1p(exp(-temp));
}

static void log_normalize(gsl_vector* x)
{
    double sum = -DBL_MAX;
    size_t i;
    const double max_val = gsl_vector_max(x);
    gsl_vector_add_constant(x, -max_val);

    for(i = 0; i < x->size; i++)
        sum = log_sum(sum, gsl_vector_get(x, i));
    gsl_vector_add_constant(x, -sum);
}

double lcfit_bsm_log_like(const double t, const bsm_t* m)
{
    if (t == 0.0 && m->b == 0.0 && m->c > m->m) {
        return -INFINITY;
    }
    else if (t == INFINITY) {
        return log(0.5) * (m->c + m->m);
    }

    double expterm = exp(-m->r * (t + m->b));
    return (m->c * log((1 + expterm) / 2) + m->m * log((1 - expterm) / 2));
}

/* The ML branch length for c, m, r, b */
double lcfit_bsm_ml_t(const bsm_t* m)
{
    if (m->c <= m->m) {
        return INFINITY;
    }

    double t = ((log((m->c - m->m) / (m->c + m->m))) / (-m->r)) - m->b;
    return t < 0.0 ? 0.0 : t;
}

double lcfit_bsm_infl_t(const bsm_t* m)
{
    lcfit_regime regime = lcfit_bsm_regime(m);

    if (!(regime == LCFIT_REGIME_1 || regime == LCFIT_REGIME_2)) {
        return NAN;
    }

    double t = -(m->b) + (1.0 / m->r) *
               log(pow(sqrt(m->c) + sqrt(m->m), 2.0) / (m->c - m->m));

    return t;
}

lcfit_regime lcfit_bsm_regime(const bsm_t* m)
{
    if (m->c == m->m) { return LCFIT_REGIME_UNKNOWN; }
    if (m->c < m->m) { return LCFIT_REGIME_4; }

    double lhs = exp(m->b * m->r);
    double rhs = pow(sqrt(m->c) + sqrt(m->m), 2.0) / (m->c - m->m);

    if (lhs > rhs) { return LCFIT_REGIME_3; }
    if (m->b > 0.0) { return LCFIT_REGIME_2; }

    return LCFIT_REGIME_1;
}

void lcfit_bsm_gradient(const double t, const bsm_t* m, double* grad)
{
    const double u = exp(-m->r * (t + m->b));

    grad[0] = log((1 + u) / 2); /* df/dc */
    grad[1] = log((1 - u) / 2); /* df/dm */
    grad[2] = (t + m->b) * (-m->c * u / (1 + u) + m->m * u / (1 - u)); /* df/dr */
    grad[3] = m->r * (-m->c * u / (1 + u) + m->m * u / (1 - u)); /* df/db */
}

/*
 * The scaling parameter for c and m to obtain log-likelihood value `l` at branch length `t`,
 * keeping `r` and `b` fixed.
 */
double lcfit_bsm_scale_factor(const double t, const double l, const bsm_t* m)
{
    double result = l / lcfit_bsm_log_like(t, m);
    return result;
}

void lcfit_bsm_rescale(const double t, const double l, bsm_t* m)
{
    double scale = lcfit_bsm_scale_factor(t, l, m);
    m->c *= scale;
    m->m *= scale;
}

/* Next, the data to fit such a log likelihood function. */
struct data_to_fit {
    const size_t n;   /* Number of observations */
    const double* t;  /* Branch lengths */
    const double* l;  /* Corresponding likelihoods */
    const double* w;  /* Corresponding weights */
    size_t iterations;
};


int lcfit_weights(const void* data, gsl_vector* weight)
{
    const size_t n = ((const struct data_to_fit*) data)->n;
    const double* l = ((const struct data_to_fit*) data)->l;
    size_t i;

    for(i = 0; i < n; i++) {
        gsl_vector_set(weight, i, l[i]);
    }

    gsl_vector_scale(weight, 1 / LAMBDA);  /* multiplies the elements of vector a by the constant factor x. */
    log_normalize(weight);

    /* for(i = 0; i < n; i++) */
    /*     fprintf(stderr, "%.3f\t", exp(gsl_vector_get(weight, i))); */
    /* fprintf(stderr, "\n"); */

    return GSL_SUCCESS;
}

void print_state_gsl(size_t iter, gsl_multifit_fdfsolver* s)
{
    gsl_vector* grad = gsl_vector_alloc(4);
    gsl_multifit_gradient(s->J, s->f, grad);

    fprintf(stderr, "G[%4zu] rsse = %.3f", iter, gsl_blas_dnrm2(s->f));
    fprintf(stderr, ", model = { %.3f, %.3f, %.6f, %.6f }",
            gsl_vector_get(s->x, 0),
            gsl_vector_get(s->x, 1),
            gsl_vector_get(s->x, 2),
            gsl_vector_get(s->x, 3));
    fprintf(stderr, ", grad = { %.6f, %.6f, %.6f, %.6f }",
            gsl_vector_get(grad, 0),
            gsl_vector_get(grad, 1),
            gsl_vector_get(grad, 2),
            gsl_vector_get(grad, 3));
    fprintf(stderr, "\n");

    gsl_vector_free(grad);
}

/* Evaluate the likelihood curve described in data at the point x. */
int lcfit_pair_f(const gsl_vector* x, void* data, gsl_vector* f)
{
    const size_t n = ((struct data_to_fit*) data)->n;
    const double* t = ((struct data_to_fit*) data)->t;
    const double* l = ((struct data_to_fit*) data)->l;
    const double* w = ((struct data_to_fit*) data)->w;
    size_t i;
    bsm_t m = {gsl_vector_get(x, 0),
               gsl_vector_get(x, 1),
               gsl_vector_get(x, 2),
               gsl_vector_get(x, 3)};

    for(i = 0; i < n; i++) {
        const double err = lcfit_bsm_log_like(t[i], &m) - l[i];

        gsl_vector_set(f, i, w[i] * err);
        /*gsl_vector_set(f, i, err * 0.25);*/
    }

    return GSL_SUCCESS;
}

/* The corresponding Jacobian. */
int lcfit_pair_df(const gsl_vector* x, void* data, gsl_matrix* J)
{
    const size_t n = ((struct data_to_fit*) data)->n;
    const double* t = ((struct data_to_fit*) data)->t;
    const double* w = ((struct data_to_fit*) data)->w;

    /* double *l = ((struct data_to_fit *) data)->l; */

    double c = gsl_vector_get(x, 0);
    double m = gsl_vector_get(x, 1);
    double r = gsl_vector_get(x, 2);
    double b = gsl_vector_get(x, 3);

    bsm_t model = {c, m, r, b};
    double grad_i[4];

    for (size_t i = 0; i < n; i++) {
        /* nx4 Jacobian matrix J(i,j) = dfi / dxj, */
        /* where fi = c*log((1+exp(-r*t[i]))/2)+m*log((1-exp(-r*t[i]))/2) - l[i] */
        /* so df/dc = log((1+exp(-r*(t+b)))/2) */
        /* so df/dm = log((1-exp(-r*(t+b)))/2) */
        /* so df/dr = c*(-(t+b))*exp(-r*(t+b))/(1+exp(-r*(t+b)))+m*(t+b)*exp(-r*(t+b))/(1-exp(-r*(t+b))) */
        /* so df/db = c*(-r)*exp(-r*(t+b))/(1+exp(-r*(t+b)))+m*r*exp(-r*(t+b))/(1-exp(-r*(t+b))) */
        /* and the xj are the parameters (c, m, r, b) */

        lcfit_bsm_gradient(t[i], &model, grad_i);

        gsl_matrix_set(J, i, 0, w[i] * grad_i[0]); /* df/dc */
        gsl_matrix_set(J, i, 1, w[i] * grad_i[1]); /* df/dm */
        gsl_matrix_set(J, i, 2, w[i] * grad_i[2]); /* df/dr */
        gsl_matrix_set(J, i, 3, w[i] * grad_i[3]); /* df/db */
    }

    return GSL_SUCCESS;
}

int lcfit_pair_fdf(const gsl_vector* x, void* data, gsl_vector* f, gsl_matrix* J)
{
    lcfit_pair_f(x, data, f);
    lcfit_pair_df(x, data, J);

    return GSL_SUCCESS;
}

void print_state_nlopt(size_t iter,
                       double sum_sq_err,
                       const double* x,
                       const double* grad)
{
    fprintf(stderr, "N[%4zu] rsse = %.3f", iter, sqrt(sum_sq_err));
    fprintf(stderr, ", model = { %.3f, %.3f, %.6f, %.6f }",
            x[0], x[1], x[2], x[3]);
    if (grad) {
        fprintf(stderr, ", grad = { %.6f, %.6f, %.6f, %.6f }",
                grad[0], grad[1], grad[2], grad[3]);
    }
    fprintf(stderr, "\n");
}

/** Least-squares objective function for fitting with NLopt.
 *
 * \param[in]  p     Number of model parameters.
 * \param[in]  x     Model parameters to evaluate.
 * \param[out] grad  Gradient of the objective function at \c x.
 * \param[in]  data  Observed likelihood data to fit.
 *
 * \return Sum of squared error from observed likelihoods
 */
double bsm_fit_objective(unsigned p,
                         const double* x,
                         double* grad,
                         void* data)
{
    struct data_to_fit* fit_data = (struct data_to_fit*) data;

    const size_t n = fit_data->n;
    const double* t = fit_data->t;
    const double* l = fit_data->l;
    const double* w = fit_data->w;

    const double c = x[0];
    const double m = x[1];
    const double r = x[2];
    const double b = x[3];

    double sum_sq_err = 0.0;

    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.0;
        grad[2] = 0.0;
        grad[3] = 0.0;
    }

    bsm_t model = {c, m, r, b};
    double grad_i[4];

    for (size_t i = 0; i < n; ++i) {
        const double err = l[i] - lcfit_bsm_log_like(t[i], &model);

        sum_sq_err += w[i] * pow(err, 2.0);

        if (grad) {
            lcfit_bsm_gradient(t[i], &model, grad_i);

            grad[0] -= 2 * w[i] * err * grad_i[0];
            grad[1] -= 2 * w[i] * err * grad_i[1];
            grad[2] -= 2 * w[i] * err * grad_i[2];
            grad[3] -= 2 * w[i] * err * grad_i[3];
        }
    }

#ifdef LCFIT4_VERBOSE
    print_state_nlopt(fit_data->iterations, sum_sq_err, x, grad);
#endif /* LCFIT4_VERBOSE */
    ++fit_data->iterations;
    return sum_sq_err;
}

int lcfit_fit_bsm(const size_t n, const double* t, const double* l, bsm_t *m, int max_iter)
{
    double *w = calloc(n, sizeof(double));
    int status, i;

    for (i = 0; i < n; i++)
	w[i] = 1.0L;
    status = lcfit_fit_bsm_weight(n, t, l, w, m, max_iter);
    free(w);
    return(status);
}

// Declare our implementations before the delegator function definition.
int lcfit_fit_bsm_weighted_gsl(const size_t, const double*, const double*, const double*, bsm_t*, int);
int lcfit_fit_bsm_weighted_nlopt(const size_t, const double*, const double*, const double*, bsm_t*, int);

int check_model(const bsm_t* m)
{
    if (m->c < 1.0 || m->m < 1.0 || m->c < m->m || m->r <= 0.0 || m->b < 0.0) {
        return 1;
    }

    if (m->r > BSM_R_MAX) {
        return 2;
    }

    return 0;
}

/**
 * The fitting procedure first tries unconstrained Levenberg-Marquardt
 * optimization to fit the model parameters to the data. See <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Nonlinear-Least_002dSquares-Fitting.html">
 * GSL's documentation</a> for method details.
 *
 * If the Levenberg-Marquardt algorithm fails to converge, returns an
 * error, or returns an invalid model, the function tries constrained
 * optimization using the SLSQP algorithm. See <a
 * href="http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#SLSQP">NLopt's
 * documentation</a> for method details.
 *
 * SLSQP is started differently depending on the reason the LM
 * algorithm failed. If LM returned a valid model, but did not
 * indicate success (e.g., failure to converge), SLSQP is used to
 * refine that model. On the other hand, if LM returned an invalid
 * model, or indicated complete failure, SLSQP starts over with the
 * original initial conditions.
 *
 * A model is considered invalid if any of the following are true:
 * - <c>m.c < 1</c>
 * - <c>m.m < 1</c>
 * - <c>m.c < m.m</c> (regime 4; see below)
 * - <c>m.r <= 0</c> (we assumed that the mutation rate is greater than zero)
 * - <c>m.b < 0</c> (branch length offset must be non-negative)
 *
 * Regime 4 is considered invalid in this scenario because, assuming a
 * tree with finite branch lengths, the probability of having more
 * mutated sites than constant sites approaches zero as sequences
 * become long.
 *
 */
int lcfit_fit_bsm_weight(const size_t n,
                         const double* t,
                         const double* l,
                         const double *w,
                         bsm_t *m,
                         int max_iter)
{
    if (n < 4) {
        fprintf(stderr, "ERROR: fitting a model requires at least four points\n");
        return LCFIT_ERROR;
    }

    bsm_t initial_model = *m;

    int status = lcfit_fit_bsm_weighted_gsl(n, t, l, w, m, max_iter);
    if (check_model(m) != 0) {
        /* GSL returned a bad model, so start over. */
        *m = initial_model;
        status = lcfit_fit_bsm_weighted_nlopt(n, t, l, w, m, max_iter);
    } else if (status != LCFIT_SUCCESS) {
        /* GSL returned a valid model but did not indicate success, so
         * try and refine the model with NLopt. */
        status = lcfit_fit_bsm_weighted_nlopt(n, t, l, w, m, max_iter);
    }

    return status;
}

int lcfit_fit_bsm_weighted_gsl(const size_t n,
                               const double* t,
                               const double* l,
                               const double *w,
                               bsm_t *m,
                               int max_iter)
{
    double x[4] = {m->c, m->m, m->r, m->b};
    int status = GSL_SUCCESS;

    struct data_to_fit d = { n, t, l, w, 0 };
    gsl_multifit_function_fdf fdf;

    /* Storing the contents of x on the stack.
     * http://www.gnu.org/software/gsl/manual/html_node/Vector-views.html */
    gsl_vector_const_view x_view = gsl_vector_const_view_array(x, 4);

    fdf.f = &lcfit_pair_f;
    fdf.df = &lcfit_pair_df;
    fdf.fdf = &lcfit_pair_fdf;
    fdf.n = n;
    fdf.p = 4; /* 4 parameters */
    fdf.params = &d;

    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver* s = gsl_multifit_fdfsolver_alloc(T, n, 4);
    assert(s != NULL && "Solver allocation failed!");
    gsl_multifit_fdfsolver_set(s, &fdf, &x_view.vector); /* Taking address of view.vector gives a const gsl_vector * */

#ifdef LCFIT4_VERBOSE
    print_state_gsl(0, s);
#endif /* LCFIT4_VERBOSE */

    do {
        d.iterations++;
        status = gsl_multifit_fdfsolver_iterate(s);

#ifdef LCFIT4_VERBOSE
        print_state_gsl(d.iterations, s);
#endif /* LCFIT4_VERBOSE */

        if (status) {
            break;
        }
        status = gsl_multifit_test_delta(s->dx, s->x, 0.0, 1e-4);
    } while (status == GSL_CONTINUE && d.iterations < max_iter);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
#ifdef LCFIT4_NOISY
    gsl_matrix* covar = gsl_matrix_alloc(4, 4);
    gsl_multifit_covar(s->J, 0.0, covar);
    gsl_matrix_fprintf(stderr, covar, "%g");

    fprintf(stderr, "c = %.5f +/- %.5f\n", FIT(0), ERR(0));
    fprintf(stderr, "m = %.5f +/- %.5f\n", FIT(1), ERR(1));
    fprintf(stderr, "r = %.5f +/- %.5f\n", FIT(2), ERR(2));
    fprintf(stderr, "b = %.5f +/- %.5f\n", FIT(3), ERR(3));

    gsl_matrix_free(covar);
#endif /* LCFIT4_NOISY */

#ifdef LCFIT4_VERBOSE
    fprintf(stderr, "[G] status = %s (%d)   iterations %zu\n",
            gsl_strerror(status), status, d.iterations);
#endif /* LCFIT4_VERBOSE */

    // translate from GSL status to LCFIT status
    // GSL error codes are defined in gsl_errno.h
    // corresonding lcfit error codes can be found in lcfit.h
    if (d.iterations >= max_iter)
        status = LCFIT_MAXITER;
    else if (status == GSL_SUCCESS)
        status = LCFIT_SUCCESS;
    else if (status ==  GSL_ENOPROG)
        status = LCFIT_ENOPROG;
    else if (status == GSL_ETOLF)
        status = LCFIT_ETOLF;
    else if (status == GSL_ETOLG)
        status = LCFIT_ETOLG;
    else
        status = LCFIT_ERROR;

    // Update fit
    m->c = FIT(0);
    m->m = FIT(1);
    m->r = FIT(2);
    m->b = FIT(3);
#undef FIT
#undef ERR

    gsl_multifit_fdfsolver_free(s);
    return status;
}

const char* nlopt_strerror(int status)
{
    switch (status) {
    case NLOPT_SUCCESS:
        return "success";
    case NLOPT_STOPVAL_REACHED:
        return "stopval reached";
    case NLOPT_FTOL_REACHED:
        return "ftol reached";
    case NLOPT_XTOL_REACHED:
        return "xtol reached";
    case NLOPT_MAXTIME_REACHED:
        return "maxtime reached";
    case NLOPT_MAXEVAL_REACHED:
        return "maxeval reached";
    case NLOPT_FAILURE:
        return "failure";
    case NLOPT_INVALID_ARGS:
        return "invalid args";
    case NLOPT_OUT_OF_MEMORY:
        return "out of memory";
    case NLOPT_ROUNDOFF_LIMITED:
        return "roundoff limited";
    case NLOPT_FORCED_STOP:
        return "forced stop";
    default:
        return "unknown status code";
    }
}

int lcfit_fit_bsm_weighted_nlopt(const size_t n,
                                 const double* t,
                                 const double* l,
                                 const double *w,
                                 bsm_t *m,
                                 int max_iter)
{
    struct data_to_fit fit_data = { n, t, l, w, 0 };

    const double lower_bounds[4] = { 1, 1, BSM_R_MIN, BSM_B_MIN };
    const double upper_bounds[4] = { INFINITY, INFINITY, BSM_R_MAX, INFINITY };

    nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, 4);
    nlopt_set_min_objective(opt, bsm_fit_objective, &fit_data);
    nlopt_set_lower_bounds(opt, lower_bounds);
    nlopt_set_upper_bounds(opt, upper_bounds);

    nlopt_set_xtol_rel(opt, 1e-4);
    nlopt_set_maxeval(opt, max_iter);

    double x[4] = { m->c, m->m, m->r, m->b };
    double minf = 0.0;

    int status = nlopt_optimize(opt, x, &minf);

#ifdef LCFIT4_VERBOSE
    fprintf(stderr, "[N] status = %s (%d)   iterations %zu\n",
            nlopt_strerror(status), status, fit_data.iterations);
#endif /* LCFIT4_VERBOSE */

    switch (status) {
    case NLOPT_SUCCESS:
    case NLOPT_STOPVAL_REACHED:
    case NLOPT_FTOL_REACHED:
    case NLOPT_XTOL_REACHED:
    case NLOPT_MAXTIME_REACHED:
        status = LCFIT_SUCCESS;
        break;
    case NLOPT_MAXEVAL_REACHED:
        status = LCFIT_MAXITER;
        break;
    case NLOPT_FAILURE:
    case NLOPT_INVALID_ARGS:
    case NLOPT_OUT_OF_MEMORY:
    case NLOPT_ROUNDOFF_LIMITED:
    case NLOPT_FORCED_STOP:
    default:
        status = LCFIT_ERROR;
    }

    m->c = x[0];
    m->m = x[1];
    m->r = x[2];
    m->b = x[3];

    nlopt_destroy(opt);
    return status;
}

#ifdef NOTYET
/* Here we minimize the KL divergence, rather than nonlinear least squares */
double kl_divergence(const double* unnorm_log_p1,
                     const double* unnorm_log_p2,
                     const size_t n)
{
    assert(unnorm_log_p1 != NULL && "Null v1");
    assert(unnorm_log_p2 != NULL && "Null v2");
    size_t i;
    gsl_vector_const_view v1 = gsl_vector_const_view_array(unnorm_log_p1, n);
    gsl_vector_const_view v2 = gsl_vector_const_view_array(unnorm_log_p2, n);

    gsl_vector* p1 = gsl_vector_alloc(n);
    gsl_vector_memcpy(p1, &v1.vector);
    log_normalize(p1);

    gsl_vector* p2 = gsl_vector_alloc(n);
    gsl_vector_memcpy(p2, &v2.vector);
    log_normalize(p2);

    gsl_vector* lr = gsl_vector_alloc(n);
    gsl_vector_memcpy(lr, p1);
    gsl_vector_sub(lr, p2);
    for(i = 0; i < n; i++) {
        gsl_vector_set(p1, i, exp(gsl_vector_get(p1, i)));
    }
    gsl_vector_mul(p1, lr);

    double kl = 0.0;
    for(i = 0; i < n; i++) {
        kl += gsl_vector_get(p1, i);
    }

    gsl_vector_free(p1);
    gsl_vector_free(p2);
    gsl_vector_free(lr);

    return kl / log(2.0);
}

static double kl_divergence_f(unsigned n, const double *x, double* grad, void *data)
{
    const double* t = ((struct data_to_fit*) data)->t;
    const double* l = ((struct data_to_fit*) data)->l;
    size_t i;

    bsm_t m = {x[0], x[1], x[2], x[3]};

    double *fit = malloc(sizeof(double) * n);
    for(i = 0; i < n; i++) {
        fit[i] = lcfit_bsm_log_like(t[i], &m);
    }

    const double kl = kl_divergence(l, fit, n);
    if(isnan(kl))
        fprintf(stderr, "NaN KL for c=%f m=%f r=%f b=%f\n", m.c, m.m, m.r, m.b);

    free(fit);

    /*fprintf(stderr, "%f\n", kl);*/

    return kl;
}

/*
 *void print_norm(const size_t n, const double* d)
 *{
 *    gsl_vector *x = gsl_vector_alloc(n);
 *    size_t i;
 *    for(i = 0; i < n; i++)
 *        gsl_vector_set(x, i, d[i]);
 *    log_normalize(x);
 *    fprintf(stderr, "%f", gsl_vector_get(x, 0));
 *    for(i = 1; i < n; i++)
 *        fprintf(stderr, "\t%f", gsl_vector_get(x, i));
 *}
 */

int lcfit_bsm_minimize_kl(const size_t n, const double* t, const double* l, bsm_t *m)
{
    /*lcfit_fit_bsm(n, t, l, m);*/
    double x[4] = {m->c, m->m, m->r, m->b};
    const size_t p = 4; /* 4 parameters */
    /*size_t i;*/

    struct data_to_fit d = {n, t, l};

    nlopt_opt opt = nlopt_create(NLOPT_LN_BOBYQA, p);
    nlopt_set_min_objective(opt, kl_divergence_f, &d);

    nlopt_set_xtol_rel(opt, 1e-2);
    nlopt_set_ftol_rel(opt, 1e-2);
    double lb[4] = {1, 1, 1e-7, 1e-7};
    double ub[4] = {DBL_MAX, DBL_MAX, DBL_MAX, 10};
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    double minf;

    int res = nlopt_optimize(opt, x, &minf);
    m->c = x[0];
    m->m = x[1];
    m->r = x[2];
    m->b = x[3];

    /*double *fit = malloc(sizeof(double) * n);*/
    /*for(i = 0; i < n; i++)*/
        /*fit[i] = lcfit_bsm_log_like(t[i], m);*/
    /*fprintf(stderr, "True: ");*/
    /*print_norm(n, l);*/
    /*fprintf(stderr, "\nFit: ");*/
    /*print_norm(n, fit);*/
    /*fprintf(stderr, "\n%f\n\n", minf);*/
    /*free(fit);*/

    nlopt_destroy(opt);

    return res < 0 ? res : 0;
}
#endif /* NOTYET */

typedef struct fn_wrapper {
    double (*fn)(double, void*);
    void* fn_args;
} fn_wrapper_t;

static double invert_wrapped_fn(double t, void* data)
{
    fn_wrapper_t* wrapper = (fn_wrapper_t*) data;
    return -(wrapper->fn(t, wrapper->fn_args));
}

// This function attempts to bisect the range [min_t, max_t] of a
// supplied function until the evaluated points enclose a maximum. If
// successful, the function updates min_t and max_t and returns true.
//
// It is assumed that the supplied function behaves similarly to one
// of the lcfit regimes. It's not guaranteed to do anything useful
// otherwise.
//
// If a minimum is encountered instead of a maximum, the function
// aborts and returns false. The function also returns false if the
// number of bisection iterations exceeds a hard-coded threshold. In
// any case, min_t and max_t are not updated if false is returned.
//
// TODO: [efficiency] likelihood function maxima are far more likely
// to be found close to zero. it would probably be more efficient if
// the guesses were biased to the left.
//
// TODO: [efficiency] the final values of t and f(t) could be passed
// back to the caller for reuse in initializing the minimizer (see
// gsl_min_fminimizer_set_with_values).

static bool bracket_maximum(double (*fn)(double, void*), void* fn_args,
                            double* min_t, double* max_t)
{
    double t[3] = {*min_t, (*min_t + *max_t) / 2.0, *max_t};
    double f[3] = {fn(t[0], fn_args), fn(t[1], fn_args), fn(t[2], fn_args)};

    const size_t MAX_ITER = 30;
    size_t iter = 0;
    bool success = false;

    for (; iter < MAX_ITER; ++iter) {
        if (f[1] > f[0] && f[1] > f[2]) {  // maximum enclosed
            success = true;
            break;
        } else if (f[1] < f[0] && f[1] < f[2]) {  // minimum enclosed
            // success = false;
            break;
        } else if (f[0] < f[1] && f[1] < f[2]) {  // monotonically increasing
            t[0] = t[1];
            f[0] = f[1];
        } else if ((f[0] > f[1] && f[1] > f[2])  // monotonically decreasing
                   || (f[1] == f[2])) {          // asymptote
            // If the rightmost points are equal, we assume they're
            // out on the asymptote to within the precision of a
            // double. Since we're also assuming the function behaves
            // similarly to one of the lcfit regimes, we won't find a
            // local maximum out there, so we bisect to the left
            // instead.

            t[2] = t[1];
            f[2] = f[1];
        }

        t[1] = (t[0] + t[2]) / 2.0;
        f[1] = fn(t[1], fn_args);
    }

#ifdef LCFIT_AUTO_VERBOSE
    fprintf(stderr, "bracket_maximum: %zu iterations\n", iter);
#endif /* LCFIT_AUTO_VERBOSE */

    if (success) {
        *min_t = t[0];
        *max_t = t[2];
    }

    return success;
}

// This function estimates the first and second derivatives of a
// well-behaved function at a point x using a five-point stencil.
//
// TODO: [efficiency] if f(x) is known beforehand, one function
// evaluation could be saved by passing that value in instead of
// recomputing it.

static void estimate_derivatives(double (*fn)(double, void*), void* fn_args,
                                 double x, double* d1, double* d2)
{
    // the central differences below are fourth order, so use a step
    // size relative to the fourth root of DBL_EPSILON
    const double h = x * pow(DBL_EPSILON, 0.25);

    // https://en.wikipedia.org/wiki/Five-point_stencil
    // https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Tables_of_selected_convolution_coefficients

    double fm2 = fn(x - 2*h, fn_args);
    double fm1 = fn(x - h, fn_args);
    double f0 = fn(x, fn_args);
    double fp1 = fn(x + h, fn_args);
    double fp2 = fn(x + 2*h, fn_args);

    *d1 = (-fp2 + 8*fp1 - 8*fm1 + fm2) / (12*h);
    *d2 = (-fp2 + 16*fp1 - 30*f0 + 16*fm1 - fm2) / (12*h*h);
}

static double find_maximum(double (*fn)(double, void*), void* fn_args,
                           double guess, double min_t, double max_t)
{
#ifdef LCFIT_AUTO_VERBOSE
    fprintf(stderr, "min = %g, guess = %g, max = %g\n", min_t, guess, max_t);
    fprintf(stderr, "f(min) = %g, f(guess) = %g, f(max) = %g\n",
            fn(min_t, fn_args),
            fn(guess, fn_args),
            fn(max_t, fn_args));
#endif /* LCFIT_AUTO_VERBOSE */

    fn_wrapper_t wrapper;
    wrapper.fn = fn;
    wrapper.fn_args = fn_args;

    gsl_function F;
    F.function = &invert_wrapped_fn;
    F.params = &wrapper;

    gsl_min_fminimizer* s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    gsl_min_fminimizer_set(s, &F, guess, min_t, max_t);

    const int MAX_ITER = 100;
    int iter = 0;
    int status;

    do {
        gsl_min_fminimizer_iterate(s);
        ++iter;

        guess = gsl_min_fminimizer_x_minimum(s);
        min_t = gsl_min_fminimizer_x_lower(s);
        max_t = gsl_min_fminimizer_x_upper(s);

        status = gsl_min_test_interval(min_t, max_t, 0.0, pow(DBL_EPSILON, 0.25));
    } while (status == GSL_CONTINUE && iter < MAX_ITER);

#ifdef LCFIT_AUTO_VERBOSE
    fprintf(stderr, "lcfit_maximize: %d iterations\n", iter);
#endif /* LCFIT_AUTO_VERBOSE */

    if (iter == MAX_ITER) {
        fprintf(stderr, "WARNING: maximum number of iterations reached during minimization\n");
    }

    gsl_min_fminimizer_free(s);

    return guess;
}

double lcfit_maximize(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                      double min_t, double max_t, double* d1, double* d2)
{
    bool is_bracketed = bracket_maximum(lnl_fn, lnl_fn_args, &min_t, &max_t);
    double guess = (min_t + max_t) / 2.0;

    if (is_bracketed) {
        guess = find_maximum(lnl_fn, lnl_fn_args, guess, min_t, max_t);
    }

    if (d1 && d2) {
        estimate_derivatives(lnl_fn, lnl_fn_args, guess, d1, d2);
    }

    return guess;
}
