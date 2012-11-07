#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

struct data {
  size_t n;
  double * y;
  double * sigma;
};

int
expb_f (const gsl_vector * x, void *params,
        gsl_vector * f)
{
  size_t n = ((struct data *)params)->n;
  double *y = ((struct data *)params)->y;
  double *sigma = ((struct data *) params)->sigma;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);
  double b = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Model Yi = A * exp(-lambda * i) + b */
      double t = i;
      double Yi = A * exp (-lambda * t) + b;
      gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
    }

  return GSL_SUCCESS;
}

int
expb_df (const gsl_vector * x, void *params,
         gsl_matrix * J)
{
  size_t n = ((struct data *)params)->n;
  double *sigma = ((struct data *) params)->sigma;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      double t = i;
      double s = sigma[i];
      double e = exp(-lambda * t);
      gsl_matrix_set (J, i, 0, e/s);
      gsl_matrix_set (J, i, 1, -t * A * e/s);
      gsl_matrix_set (J, i, 2, 1/s);

    }
  return GSL_SUCCESS;
}

int
expb_fdf (const gsl_vector * x, void *params,
          gsl_vector * f, gsl_matrix * J)
{
  expb_f (x, params, f);
  expb_df (x, params, J);

  return GSL_SUCCESS;
}

//The main part of the program sets up a Levenberg-Marquardt solver and some simulated random data. The data uses the known parameters (1.0,5.0,0.1) combined with gaussian noise (standard deviation = 0.1) over a range of 40 timesteps. The initial guess for the parameters is chosen as (0.0, 1.0, 0.0).

int
main (void)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 40;
  const size_t p = 3;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);

  double y[n], sigma[n];

  struct data d = { n, y, sigma};

  gsl_multifit_function_fdf f;

  double x_init[3] = { 1.0, 0.0, 0.0 };

  gsl_vector_view x = gsl_vector_view_array (x_init, p);

  const gsl_rng_type * type;
  gsl_rng * r;

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  /* This is the data to be fitted */

  for (i = 0; i < n; i++)
    {
      double t = i;
      y[i] = 1.0 + 5 * exp (-0.1 * t)
                 + gsl_ran_gaussian(r, 0.1);
      sigma[i] = 0.1;
      printf("data: %d %g %g\n", i, y[i], sigma[i]);
    };

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      printf ("status = %s\n", gsl_strerror (status));

      print_state (iter, s);

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);

  gsl_multifit_covar (s->J, 0.0, covar);

  gsl_matrix_fprintf (stdout, covar, "%g");

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  printf("A      = %.5f +/- %.5f\n", FIT(0), ERR(0));
  printf("lambda = %.5f +/- %.5f\n", FIT(1), ERR(1));
  printf("b      = %.5f +/- %.5f\n", FIT(2), ERR(2));

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multifit_fdfsolver_free (s);
  return 0;
}

int
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_blas_dnrm2 (s->f));
}

//The iteration terminates when the change in x is smaller than 0.0001, as both an absolute and relative change. Here are the results of running the program,
//
//iter: 0 x = 1.00000000 0.00000000 0.00000000 |f(x)| = 118.574
//iter: 1 x = 1.64919392 0.01780040 0.64919392 |f(x)| = 77.2068
//iter: 2 x = 2.86269020 0.08032198 1.45913464 |f(x)| = 38.0579
//iter: 3 x = 4.97908864 0.11510525 1.06649948 |f(x)| = 10.1548
//iter: 4 x = 5.03295496 0.09912462 1.00939075 |f(x)| = 6.4982
//iter: 5 x = 5.05811477 0.10055914 0.99819876 |f(x)| = 6.33121
//iter: 6 x = 5.05827645 0.10051697 0.99756444 |f(x)| = 6.33119
//iter: 7 x = 5.05828006 0.10051819 0.99757710 |f(x)| = 6.33119
//
//A      = 5.05828 +/- 0.05983
//lambda = 0.10052 +/- 0.00309
//b      = 0.99758 +/- 0.03944
//status = success
//
//The approximate values of the parameters are found correctly. The errors on the parameters are given by the square roots of the diagonal elements of the covariance matrix.
