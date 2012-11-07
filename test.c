#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

struct data {
    size_t n;
    double * t;
    double * l;
};

int expb_f(const gsl_vector * x, void *params, gsl_vector * f)
{
    size_t n = ((struct data *)params)->n;
    double *t = ((struct data *)params)->t;
    double *l = ((struct data *) params)->l;

    double c = gsl_vector_get(x, 0);
    double m = gsl_vector_get(x, 1);
    double r = gsl_vector_get(x, 2);

    size_t i;
    double ert;

    for(i = 0; i < n; i++) {
        /* Model l[i] = c*log((1+exp(-r*t[i]))/2)+m*log((1-exp(-r*t[i]))/2) */
        ert = exp(-r * t[i]);
        gsl_vector_set(f, i, c * log((1 + ert) / 2) + m * log((1 - ert) / 2) - l[i]);
    }

    return GSL_SUCCESS;
}

int expb_df(const gsl_vector * x, void *params, gsl_matrix * J)
{
    size_t n = ((struct data *)params)->n;
    double *t = ((struct data *)params)->t;
    double *l = ((struct data *) params)->l;

    double c = gsl_vector_get(x, 0);
    double m = gsl_vector_get(x, 1);
    double r = gsl_vector_get(x, 2);

    size_t i;
    double ert;

    for(i = 0; i < n; i++) {
        /* nx3 Jacobian matrix J(i,j) = dfi / dxj, */
        /* where fi = = c*log((1+exp(-r*t[i]))/2)+m*log((1-exp(-r*t[i]))/2) - l[i] */
        /* and the xj are the parameters (c, m, r) */

        ert = exp(-r * t[i]);
        gsl_matrix_set(J, i, 0, log((1 + ert) / 2)); /* df/dc */
        gsl_matrix_set(J, i, 1, log((1 - ert) / 2)); /* df/dm */
        gsl_matrix_set(J, i, 2, t[i] * (-c * ert / (1 + ert) + m * ert / (1 - ert))); /* df/dt */
    }
    return GSL_SUCCESS;
}

int expb_fdf(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
    expb_f(x, params, f);
    expb_df(x, params, J);

    return GSL_SUCCESS;
}

int print_state(unsigned int iter, gsl_multifit_fdfsolver * s)
{
    printf("iter: %3u x = % 15.8f % 15.8f % 15.8f "
           "|f(x)| = %g\n",
           iter,
           gsl_vector_get(s->x, 0),
           gsl_vector_get(s->x, 1),
           gsl_vector_get(s->x, 2),
           gsl_blas_dnrm2(s->f));
}

int main(void)
{
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;

    int status;
    size_t i;
    unsigned int iter = 0;

    const size_t n = 3;

    gsl_matrix *covar = gsl_matrix_alloc(3, 3);

    double t[3] = {0.25, 0.5, 1}, l[3] = { -2.787877533411486, -2.721250292026327, -3.05124979115564};

    struct data d = {n, t, l};

    gsl_multifit_function_fdf f;

    double x_init[3] = { 3.0, 4.0, 5.0 };

    gsl_vector_view x = gsl_vector_view_array(x_init, 3);

    f.f = &expb_f;
    f.df = &expb_df;
    f.fdf = &expb_fdf;
    f.n = n;
    f.p = 3;
    f.params = &d;

    /* This is the data to be fitted */

    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc(T, n, 3);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);

    print_state(iter, s);

    do {
        iter++;
        status = gsl_multifit_fdfsolver_iterate(s);

        printf("status = %s\n", gsl_strerror(status));

        print_state(iter, s);

        if(status)
            break;

        status = gsl_multifit_test_delta(s->dx, s->x,
                                         1e-4, 1e-4);
    } while(status == GSL_CONTINUE && iter < 500);

    gsl_multifit_covar(s->J, 0.0, covar);

    gsl_matrix_fprintf(stdout, covar, "%g");

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    printf("c = %.5f +/- %.5f\n", FIT(0), ERR(0));
    printf("m = %.5f +/- %.5f\n", FIT(1), ERR(1));
    printf("r = %.5f +/- %.5f\n", FIT(2), ERR(2));

    printf("status = %s\n", gsl_strerror(status));

    gsl_multifit_fdfsolver_free(s);
    return 0;
}

