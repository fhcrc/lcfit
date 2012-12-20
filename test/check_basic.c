#include "lcfit.h"

#include <check.h>
#include <memory.h>
#include <math.h>
#include <gsl/gsl_errno.h>

#define TOL 1e-4

START_TEST(test_scale)
{
  const double c = 1500, m = 1000, r = 1.0, b = 0.5;
  const double t = 0.5, l = -23804.3;
  /* Try scaling */
  const double scale = lcfit_cfn_scale_factor(t, l, c, m, r, b);

  /* Test LL at t */
  const double calc_ll = lcfit_cfn_log_like(t, c*scale, m*scale, r, b);

  fail_unless(fabs(l - calc_ll) < TOL, 
              "LL after scaling was %f instead of %f", calc_ll, l);
}
END_TEST

/* Run a single fit, require that it converge, and the residuals decrase from
 * initial conditions */
void
fail_unless_fit_improves(const double x_init[4], const double t[4], const double l[4])
{
  double x[4];
  memcpy(x, x_init, 4*sizeof(double));

  int result = lcfit_fit_cfn(4, t, l, x);
  fail_unless(!result, "fit_ll returned '%s'", gsl_strerror(result));

  /* Estimates must improve */
  int i;
  double init_residuals = 0, updated_residuals = 0;
  for(i = 0; i < 4; ++i) {
    fail_unless(x[i] != x_init[i]);
    double init_fit_ll = lcfit_cfn_log_like(t[i], x_init[0], x_init[1], x_init[2], x_init[3]),
           new_fit_ll = lcfit_cfn_log_like(t[i], x[0], x[1], x[2], x[3]);
    init_residuals += fabs(init_fit_ll - l[i]);
    updated_residuals += fabs(new_fit_ll - l[i]);
  }
  fail_unless(init_residuals > updated_residuals,
      "residuals have increased after fitting (%f -> %f)",
      init_residuals, updated_residuals);
}

/* Basic test of fit_ll - just makes sure it runs */
START_TEST(test_fit_ll)
{
  /* Scaled c and m from inputs to test_scale */
  const double x_init[4] = {20739.66, 13826.44, 1.0, 0.5};
  const double t[4] = {0.1, 0.2, 0.5, 1.0};
  const double l[4] = {-23912.9, -23861.9, -23804.3, -23820.9};
  fail_unless_fit_improves(x_init, t, l);
}
END_TEST

/* Check that l == ll(t, c, m, r, b) */
void fail_unless_ll_matches(float l, float t, float c, float m, float r, float b)
{
  float res = lcfit_cfn_log_like(t, c, m, r, b);
  fail_unless(fabs(l - res) < TOL, "%f != %f", l, res);
}

/* Check a few evaluations of CFN */
START_TEST(test_cfn_ll)
{
  fail_unless_ll_matches(-8.692919, 0.2, 100, 1, 0.2, 0.4);
  fail_unless_ll_matches(-695.2381, 0.6, 1000, 250, 0.2, 0.4);
}
END_TEST

Suite*
lcfit_suite (void)
{
  Suite *s = suite_create ("LCFIT");

  /* Core test case */
  TCase *tc_core = tcase_create ("LCFIT");
  tcase_add_test (tc_core, test_scale);
  tcase_add_test (tc_core, test_cfn_ll);
  tcase_add_test (tc_core, test_fit_ll);
  suite_add_tcase (s, tc_core);

  return s;
}

int
main (void)
{
  int number_failed;
  Suite *s = lcfit_suite ();
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
