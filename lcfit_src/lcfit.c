/* http://www.gnu.org/software/gsl/manual/html_node/Nonlinear-Least_002dSquares-Fitting.html */
/**
 * \file lcfit.c
 * \brief lcfit C-API implementation.
 */

#include "lcfit.h"

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_roots.h>

#include <nlopt.h>

const static double LAMBDA = 50;

const bsm_t DEFAULT_INIT = {1500.0, 1000.0, 1.0, 0.5};

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
    double expterm = exp(-m->r * (t + m->b));
    return (m->c * log((1 + expterm) / 2) + m->m * log((1 - expterm) / 2));
}

/* The ML branch length for c, m, r, b */
double lcfit_bsm_ml_t(const bsm_t* m)
{
    double t = ((log((m->c - m->m) / (m->c + m->m))) / (-m->r)) - m->b;
    return t < 0.0 ? 0.0 : t;
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

    size_t i;
    double expterm;

    for(i = 0; i < n; i++) {
        /* nx4 Jacobian matrix J(i,j) = dfi / dxj, */
        /* where fi = c*log((1+exp(-r*t[i]))/2)+m*log((1-exp(-r*t[i]))/2) - l[i] */
        /* so df/dc = log((1+exp(-r*(t+b)))/2) */
        /* so df/dm = log((1-exp(-r*(t+b)))/2) */
        /* so df/dr = c*(-(t+b))*exp(-r*(t+b))/(1+exp(-r*(t+b)))+m*(t+b)*exp(-r*(t+b))/(1-exp(-r*(t+b))) */
        /* so df/db = c*(-r)*exp(-r*(t+b))/(1+exp(-r*(t+b)))+m*r*exp(-r*(t+b))/(1-exp(-r*(t+b))) */
        /* and the xj are the parameters (c, m, r, b) */

        expterm = exp(-r * (t[i] + b));
        gsl_matrix_set(J, i, 0, w[i] * log((1 + expterm) / 2)); /* df/dc */
        gsl_matrix_set(J, i, 1, w[i] * log((1 - expterm) / 2)); /* df/dm */
        gsl_matrix_set(J, i, 2, w[i] * ((t[i] + b) * (-c * expterm / (1 + expterm) + m * expterm / (1 - expterm)))); /* df/dr */
        gsl_matrix_set(J, i, 3, w[i] * (r * (-c * expterm / (1 + expterm) + m * expterm / (1 - expterm)))); /* df/db */
    }

    return GSL_SUCCESS;
}

int lcfit_pair_fdf(const gsl_vector* x, void* data, gsl_vector* f, gsl_matrix* J)
{
    lcfit_pair_f(x, data, f);
    lcfit_pair_df(x, data, J);

    return GSL_SUCCESS;
}

void print_state(unsigned int iter, gsl_multifit_fdfsolver* s)
{
    printf("iter: %3u {c,m,r,b} = % 15.8f % 15.8f % 15.8f % 15.8f "
           "|f(x)| = %g\n",
           iter,
           gsl_vector_get(s->x, 0),
           gsl_vector_get(s->x, 1),
           gsl_vector_get(s->x, 2),
           gsl_vector_get(s->x, 3),
           gsl_blas_dnrm2(s->f));
}


/** \brief convenience function for non-weighted lcfit
 *
 * \param n Number of observations in \c t and \c l
 * \param t Branch length
 * \param l Log-likelihood values at \c t
 * \param m Initial conditions for the model.
 * Combine #DEFAULT_INIT and #lcfit_bsm_scale_factor for reasonable starting conditions.
 *
 * @return status==0 for success, non-zero otherwise
 */
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

/** \brief Fit the BSM
 *
 * \param n Number of observations in \c t and \c l
 * \param t Branch length
 * \param l Log-likelihood values at \c t
 * \param w weight for sample point at \c t
 * \param m Initial conditions for the model.
 * Combine #DEFAULT_INIT and #lcfit_bsm_scale_factor for reasonable starting conditions.
 */
int xxx_lcfit_fit_bsm_weight(const size_t n, const double* t, const double* l, const double *w, bsm_t *m, int max_iter)
{
    double x[4] = {m->c, m->m, m->r, m->b};
    int status = GSL_SUCCESS;
    unsigned int iter = 0;

    struct data_to_fit d = {n, t, l, w};
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

    do {
        iter++;
        status = gsl_multifit_fdfsolver_iterate(s);

#ifdef VERBOSE
        printf("status = %s (%d)\n", gsl_strerror(status), status);
        print_state(iter, s);
#endif /* VERBOSE */

        if (status) {
            break;
	}
        status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
    } while(status == GSL_CONTINUE && iter < max_iter);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
#ifdef VERBOSE
    gsl_matrix* covar = gsl_matrix_alloc(4, 4);
    gsl_multifit_covar(s->J, 0.0, covar);
    gsl_matrix_fprintf(stdout, covar, "%g");

    printf("c = %.5f +/- %.5f\n", FIT(0), ERR(0));
    printf("m = %.5f +/- %.5f\n", FIT(1), ERR(1));
    printf("r = %.5f +/- %.5f\n", FIT(2), ERR(2));
    printf("b = %.5f +/- %.5f\n", FIT(3), ERR(3));

    printf("status = %s (%d)   iterations %d\n", gsl_strerror(status), status, iter);
    gsl_matrix_free(covar);
#endif /* VERBOSE */
    printf("status = %s (%d)   iterations %d\n", gsl_strerror(status), status, iter);

    // translate from GSL status to LCFIT status
    // GSL error codes are defined in gsl_errno.h
    // a local copy can be found in ~matsengrp/local/include/gsl/gsl_errno.h
    // corresonding lcfit error codes can be found in lcfit.h
    if (iter >= max_iter) {
	status = LCFIT_MAXITER;
    } else if (status == GSL_SUCCESS) 
	status = LCFIT_SUCCESS;
    else if (status ==  GSL_ENOPROG)
	status = LCFIT_ENOPROG;
    else if (status == GSL_ETOLF)
	status = LCFIT_ETOLF;
    //else if (status == GSL_ETOLG)
    //  status = LCFIT_ETOLG;
    else {
	fprintf(stderr, "GSL status - iteration: %d status: %d, %s\n", iter, status, gsl_strerror(status));
	status = LCFIT_ERROR;
    }

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

double bsm_fit_objective(unsigned p,
                         const double* x,
                         double* grad,
                         void* data)
{
    struct data_to_fit* fit_data = (struct data_to_fit*) data;

    const size_t n = fit_data->n;
    const double* t = fit_data->t;
    const double* l_hat = fit_data->l;
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

    for (size_t i = 0; i < n; ++i) {
        const double u = exp(-r * (t[i] + b));
        const double l = c * log((1 + u) / 2) +
                         m * log((1 - u) / 2);

        const double err = l_hat[i] - l;

        sum_sq_err += w[i] * pow(err, 2.0);

        if (grad) {
            grad[0] -= 2 * w[i] * err * log((1 + u) / 2);
            grad[1] -= 2 * w[i] * err * log((1 - u) / 2);
            grad[2] -= 2 * w[i] * err * (t[i] + b) * (-c * u / (1 + u)) + m * u / (1 - u);
            grad[3] -= 2 * w[i] * err * r * (-c * u / (1 + u)) + m * u / (1 - u);
        }
    }

    ++fit_data->iterations;
    return sum_sq_err;
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

int lcfit_fit_bsm_weight(const size_t n,
                         const double* t,
                         const double* l,
                         const double *w,
                         bsm_t *m,
                         int max_iter)
{
    struct data_to_fit fit_data = { n, t, l, w };

    const double lower_bounds[4] = { 1, 1, 1e-7, 1e-4 };
    const double upper_bounds[4] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, 10 };

    nlopt_opt opt = nlopt_create(NLOPT_LD_MMA, 4);
    nlopt_set_min_objective(opt, bsm_fit_objective, &fit_data);
    nlopt_set_lower_bounds(opt, lower_bounds);
    nlopt_set_upper_bounds(opt, upper_bounds);

    nlopt_set_xtol_abs1(opt, 1e-4);
    nlopt_set_xtol_rel(opt, 1e-4);
    //nlopt_set_maxeval(opt, max_iter);

    double x[4] = { m->c, m->m, m->r, m->b };
    double minf = 0.0;

    int status = nlopt_optimize(opt, x, &minf);

    switch (status) {
    case NLOPT_SUCCESS:
    case NLOPT_STOPVAL_REACHED:
    case NLOPT_FTOL_REACHED:
    case NLOPT_XTOL_REACHED:
    case NLOPT_MAXTIME_REACHED:
        status = LCFIT_SUCCESS;
        break;
    case NLOPT_MAXEVAL_REACHED:
        //printf("status = %s (%d)   iterations %zu\n",
        //       nlopt_strerror(status), status, fit_data.iterations);
        status = LCFIT_MAXITER;
        break;
    case NLOPT_FAILURE:
    case NLOPT_INVALID_ARGS:
    case NLOPT_OUT_OF_MEMORY:
    case NLOPT_ROUNDOFF_LIMITED:
    case NLOPT_FORCED_STOP:
    default:
        //printf("status = %s (%d)   iterations %zu\n",
        //       nlopt_strerror(status), status, fit_data.iterations);
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
