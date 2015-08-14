#include "catch.hpp"

#include <algorithm>
#include <iostream>
#include <gsl/gsl_rng.h>

#include "lcfit.h"
#include "lcfit_cpp.h"
#include "lcfit_rejection_sampler.h"

#define TOL 1e-4

using namespace lcfit;

TEST_CASE("retain_top", "")
{
    const std::vector<Point> points{{0, 1}, {0.1, 1.2}, {2, 0.6}};
    std::vector<Point> expected{{0.1,1.2}};
    std::vector<Point> result = retain_top(points, 1);
    REQUIRE(result.size() == 1);
    REQUIRE(result == expected);

    // Try with two points
    expected = {{0,1},{0.1,1.2}};
    result = retain_top(points, 2);
    REQUIRE(result.size() == 2);
    REQUIRE(result == expected);
}

TEST_CASE("monotonicity_decreasing", "Ensure that monotonicity correctly identifies decreasing points")
{
    REQUIRE(monotonicity({{0, 1}, {0.1, 0.8}, {2, 0.6}}) == Monotonicity::MONO_DEC);
}

TEST_CASE("monotonicity_nonmonotonic", "Ensure that monotonicity correctly points including an extremum")
{
    REQUIRE(monotonicity({{0, 1}, {0.1, 0.8}, {2, 0.81}}) == Monotonicity::NON_MONOTONIC);
}

TEST_CASE("monotonicity_inc", "Ensure that monotonicity correctly points including an extremum") {
    REQUIRE(monotonicity({{0, 1}, {0.1, 2}, {2, 3}}) == Monotonicity::MONO_INC);
}

TEST_CASE("test_scale", "")
{
  const bsm_t m = {1500, 1000, 1, 0.5};
  bsm_t scaled = m;
  REQUIRE(scaled.m == m.m);

  REQUIRE(&m != &scaled);

  const double t = 0.5, l = -23804.3;
  /* Try scaling */
  const double scale = lcfit_bsm_scale_factor(t, l, &m);

  /* Test LL at t */
  scaled.c *= scale;
  scaled.m *= scale;
  const double calc_ll = lcfit_bsm_log_like(t, &scaled);

  REQUIRE(fabs(l - calc_ll) < TOL);
}

/* Run a single fit, require that it converge, and the residuals decrase from
 * initial conditions */
void
fail_unless_fit_improves(const bsm_t* m, const double t[4], const double l[4])
{
  bsm_t fit = *m;

  int result = lcfit_fit_bsm(4, t, l, &fit, 500);
  REQUIRE(!result);

  /* Estimates must improve */
  int i;
  double init_residuals = 0, updated_residuals = 0;
  for(i = 0; i < 4; ++i) {
    double init_fit_ll = lcfit_bsm_log_like(t[i], m),
           new_fit_ll = lcfit_bsm_log_like(t[i], &fit);
    init_residuals += fabs(init_fit_ll - l[i]);
    updated_residuals += fabs(new_fit_ll - l[i]);
  }
  REQUIRE(init_residuals > updated_residuals);
}

/* Basic test of fit_ll - just makes sure it runs */
TEST_CASE("test_fit_ll", "")
{
  /* Scaled c and m from inputs to test_scale */
  const bsm_t m = {20739.66, 13826.44, 1.0, 0.5};
  const double t[4] = {0.1, 0.2, 0.5, 1.0};
  const double l[4] = {-23912.9, -23861.9, -23804.3, -23820.9};
  fail_unless_fit_improves(&m, t, l);
}

/* Check that l == ll(t, c, m, r, b) */
void
fail_unless_ll_matches(float l, float t, const bsm_t* m)
{
  float res = lcfit_bsm_log_like(t, m);
  REQUIRE(fabs(l - res) < TOL);
}

/* Check a few evaluations of BSM */
TEST_CASE("test_bsm_ll", "")
{
  bsm_t m1 = {100, 1, 0.2, 0.4};
  fail_unless_ll_matches(-8.692919, 0.2, &m1);
  bsm_t m2 = {1000, 250, 0.2, 0.4};
  fail_unless_ll_matches(-695.2381, 0.6, &m2);
}

TEST_CASE("test_bsm_fit", "Test fitting an actual BSM log-likelihood function")
{
    auto log_like = [](const double t) -> double {
        const static bsm_t m = {2000, 500, 2.0, 0.4};
        return lcfit_bsm_log_like(t, &m);
    };
    const std::vector<double> t{0.1,0.15,0.5};
    const bsm_t m = {1500, 1000, 1.0, 0.5};
    lcfit::LCFitResult r = fit_bsm_log_likelihood(log_like, m, t);
    bsm_t* fit = &r.model_fit;

    // Fit should be quite close for all evaluated points.
    for(const Point& p : r.evaluated_points) {
        REQUIRE(std::abs(lcfit_bsm_log_like(p.x, fit) - p.y) < TOL);
    }
}

TEST_CASE("test_rejection_sampler", "Test sampling from a BSM log-likelihood function")
{
    auto log_like = [](const double t) -> double {
        const static bsm_t m = {2000, 500, 2.0, 0.4};
        return lcfit_bsm_log_like(t, &m);
    };
    const std::vector<double> t{0.1,0.15,0.5};
    const bsm_t m = {1500, 1000, 1.0, 0.5};
    lcfit::LCFitResult r = fit_bsm_log_likelihood(log_like, m, t);

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

    lcfit::rejection_sampler sampler(rng, r);
    double s = sampler.sample();
    REQUIRE(s > 0.0);
}
