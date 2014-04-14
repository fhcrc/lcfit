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
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_roots.h>

const bsm_t DEFAULT_INIT = {1500.0, 1000.0, 1.0, 0.5};

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
};

/* Evaluate the likelihood curve described in data at the point x. */
int lcfit_pair_f(const gsl_vector* x, void* data, gsl_vector* f)
{
    const size_t n = ((struct data_to_fit*) data)->n;
    const double* t = ((struct data_to_fit*) data)->t;
    const double* l = ((struct data_to_fit*) data)->l;

    bsm_t m = {gsl_vector_get(x, 0),
               gsl_vector_get(x, 1),
               gsl_vector_get(x, 2),
               gsl_vector_get(x, 3)};

    size_t i;

    for(i = 0; i < n; i++) {
        gsl_vector_set(f, i, lcfit_bsm_log_like(t[i], &m) - l[i]);
    }

    return GSL_SUCCESS;
}

/* The corresponding Jacobian. */
int lcfit_pair_df(const gsl_vector* x, void* data, gsl_matrix* J)
{
    const size_t n = ((struct data_to_fit*) data)->n;
    const double* t = ((struct data_to_fit*) data)->t;
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
        gsl_matrix_set(J, i, 0, log((1 + expterm) / 2)); /* df/dc */
        gsl_matrix_set(J, i, 1, log((1 - expterm) / 2)); /* df/dm */
        gsl_matrix_set(J, i, 2, (t[i] + b) * (-c * expterm / (1 + expterm) + m * expterm / (1 - expterm))); /* df/dr */
        gsl_matrix_set(J, i, 3, r * (-c * expterm / (1 + expterm) + m * expterm / (1 - expterm))); /* df/db */
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
    printf("iter: %3u x = % 15.8f % 15.8f % 15.8f % 15.8f "
           "|f(x)| = %g\n",
           iter,
           gsl_vector_get(s->x, 0),
           gsl_vector_get(s->x, 1),
           gsl_vector_get(s->x, 2),
           gsl_vector_get(s->x, 3),
           gsl_blas_dnrm2(s->f));
}

/*
     c = x[0];
     r = x[1];
     m = x[2];
     b = x[3];
*/
int lcfit_fit_bsm(const size_t n, const double* t, const double* l, bsm_t *m)
{
    double x[4] = {m->c, m->m, m->r, m->b};

    int status;
    unsigned int iter = 0;

    struct data_to_fit d = {n, t, l};
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
        printf("status = %s\n", gsl_strerror(status));
        print_state(iter, s);
#endif /* VERBOSE */

        if(status)
            break;

        status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
    } while(status == GSL_CONTINUE && iter < 500);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
#ifdef VERBOSE
    gsl_matrix* covar = gsl_matrix_alloc(4, 4);
    gsl_multifit_covar(s->J, 0.0, covar);
    gsl_matrix_fprintf(stdout, covar, "%g");
    gsl_matrix_free(covar);

    printf("c = %.5f +/- %.5f\n", FIT(0), ERR(0));
    printf("m = %.5f +/- %.5f\n", FIT(1), ERR(1));
    printf("r = %.5f +/- %.5f\n", FIT(2), ERR(2));
    printf("r = %.5f +/- %.5f\n", FIT(3), ERR(3));

    printf("status = %s\n", gsl_strerror(status));
#endif /* VERBOSE */

    // Update fit
    m->c = FIT(0);
    m->m = FIT(1);
    m->r = FIT(2);
    m->b = FIT(3);
#undef FIT
#undef ERR

    gsl_multifit_fdfsolver_free(s);
    return 0;
}

/* Here we minimize the KL divergence, rather than nonlinear least squares */
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
    for(i = 0; i < x->size; i++)
        sum = log_sum(sum, gsl_vector_get(x, i));
    gsl_vector_add_constant(x, -sum);
}

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

static double kl_divergence_f(const gsl_vector *x, void *data)
{
    const size_t n = ((struct data_to_fit*) data)->n;
    const double* t = ((struct data_to_fit*) data)->t;
    const double* l = ((struct data_to_fit*) data)->l;
    size_t i;

    bsm_t m = {gsl_vector_get(x, 0),
               gsl_vector_get(x, 1),
               gsl_vector_get(x, 2),
               gsl_vector_get(x, 3)};

    double *fit = malloc(sizeof(double) * n);
    for(i = 0; i < n; i++) {
        fit[i] = lcfit_bsm_log_like(t[i], &m);
    }

    const double kl = kl_divergence(l, fit, n);

    free(fit);

    return kl;
}

int lcfit_bsm_minimize_kl(const size_t n, const double* t, const double* l, bsm_t *m)
{
    double x[4] = {m->c, m->m, m->r, m->b};
    const size_t p = 4; /* 4 parameters */

    int status;
    unsigned int iter = 0;

    struct data_to_fit d = {n, t, l};

    /* Storing the contents of x on the stack.
     * http://www.gnu.org/software/gsl/manual/html_node/Vector-views.html */
    gsl_vector_const_view x_view = gsl_vector_const_view_array(x, 4);

    gsl_vector *step_size = gsl_vector_alloc(4);
    gsl_vector_set_all(step_size, 0.05);
    gsl_vector_mul(step_size, &x_view.vector);

    gsl_multimin_function f;
    f.f = &kl_divergence_f;
    f.n = p;
    f.params = &d;

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *solver = gsl_multimin_fminimizer_alloc(T, p);
    assert(solver != NULL && "Solver allocation failed!");
    gsl_multimin_fminimizer_set(solver, &f, &x_view.vector, step_size);

    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(solver);

        if(status)
            break;

//#ifdef VERBOSE
        printf("status = %s\n", gsl_strerror(status));
        /*print_state(iter, s);*/
//#endif /* VERBOSE */
    } while(status == GSL_CONTINUE && iter < 500);

#define FIT(i) gsl_vector_get(solver->x, i)
//#ifdef VERBOSE

    printf("c = %.5f\n", FIT(0));
    printf("m = %.5f\n", FIT(1));
    printf("r = %.5f\n", FIT(2));
    printf("b = %.5f\n", FIT(3));

    printf("status = %s\n", gsl_strerror(status));
//#endif /* VERBOSE */

    // Update fit
    m->c = FIT(0);
    m->m = FIT(1);
    m->r = FIT(2);
    m->b = FIT(3);
#undef FIT

    gsl_multimin_fminimizer_free(solver);
    gsl_vector_free(step_size);
    return 0;
}
