#include "catch.hpp"

#include <algorithm>
#include <iostream>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>

#include "lcfit.h"
#include "lcfit_cpp.h"
#include "lcfit_rejection_sampler.h"

#define TOL 1e-4

using namespace lcfit;

const bsm_t REGIME_1 = {10.0, 1.0, 1.0, 0.0};
const bsm_t REGIME_2 = {10.0, 1.0, 1.0, 0.1};
const bsm_t REGIME_3 = {10.0, 1.0, 1.0, 1.0};
const bsm_t REGIME_4 = {1.0, 10.0, 1.0, 0.1};

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


TEST_CASE("test points are selected properly (C++)", "[select_points_cpp]") {
    bsm_t model = {1200.0, 300.0, 1.0, 0.2}; // ml_t = 0.310826
    auto log_like = [model](double t) { return lcfit_bsm_log_like(t, &model); };
    const size_t max_pts = 8;

    SECTION("when the maximum is to the left of the starting points") {
        std::vector<double> t = {0.5, 1.0, 1.1};
        std::vector<Point> selected_pts = select_points(log_like, t, max_pts);

        REQUIRE(selected_pts.size() == 4);

        REQUIRE(selected_pts[0].x == Approx(t[0] / 10.0));
        REQUIRE(selected_pts[1].x == Approx(t[0]));
        REQUIRE(selected_pts[2].x == Approx(t[1]));
        REQUIRE(selected_pts[3].x == Approx(t[2]));

        REQUIRE(selected_pts[1].y > selected_pts[0].y);
        REQUIRE(selected_pts[1].y > selected_pts[2].y);
    }

    SECTION("when the maximum is to the right of the starting points") {
        std::vector<double> t = {0.01, 0.1, 0.2};
        std::vector<Point> selected_pts = select_points(log_like, t, max_pts);

        REQUIRE(selected_pts.size() == 5);

        REQUIRE(selected_pts[0].x == Approx(t[0]));
        REQUIRE(selected_pts[1].x == Approx(t[1]));
        REQUIRE(selected_pts[2].x == Approx(t[2]));
        REQUIRE(selected_pts[3].x == Approx(2.0 * t[2]));
        REQUIRE(selected_pts[4].x == Approx(2.0 * 2.0 * t[2]));

        REQUIRE(selected_pts[3].y > selected_pts[2].y);
        REQUIRE(selected_pts[3].y > selected_pts[4].y);
    }

    SECTION("when the maximum is enclosed by the starting points") {
        std::vector<double> t = {0.1, 0.5, 1.0};
        std::vector<Point> selected_pts = select_points(log_like, t, max_pts);

        REQUIRE (selected_pts.size() == 3);

        REQUIRE(selected_pts[0].x == Approx(t[0]));
        REQUIRE(selected_pts[1].x == Approx(t[1]));
        REQUIRE(selected_pts[2].x == Approx(t[2]));

        REQUIRE(selected_pts[1].y > selected_pts[0].y);
        REQUIRE(selected_pts[1].y > selected_pts[2].y);
    }
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
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    double lambda = 0.1;

    SECTION("multiple times") {
        bsm_t m = {10.0, 1.0, 1.0, 0.0};
        lcfit::rejection_sampler sampler(rng, m, lambda);

        std::vector<double> samples = sampler.sample_n(1000);
        REQUIRE(std::all_of(samples.begin(), samples.end(),
                            [](const double x) { return x > 0.0; }));
    }

    SECTION("in regime 1") {
        lcfit::rejection_sampler sampler(rng, REGIME_1, lambda);

        double s = sampler.sample();
        REQUIRE(s > 0.0);
    }

    SECTION("in regime 2") {
        lcfit::rejection_sampler sampler(rng, REGIME_2, lambda);

        double s = sampler.sample();
        REQUIRE(s > 0.0);
    }

    SECTION("in regime 3") {
        lcfit::rejection_sampler sampler(rng, REGIME_3, lambda);

        double s = sampler.sample();
        REQUIRE(s > 0.0);
    }

    SECTION("in regime 4") {
        lcfit::rejection_sampler sampler(rng, REGIME_4, lambda);

        double s = sampler.sample();
        REQUIRE(s > 0.0);
    }

    gsl_rng_free(rng);
}

bool test_sample_distribution(const std::vector<double>& samples,
                              const lcfit::rejection_sampler& sampler,
                              const size_t n_bins, const double tolerance)
{
    auto range = std::minmax_element(samples.begin(), samples.end());

    gsl_histogram* h = gsl_histogram_alloc(n_bins);
    gsl_histogram_set_ranges_uniform(h, *(range.first), *(range.second));
    const double bin_width = (*(range.second) - *(range.first)) / n_bins;

    for (const double& s : samples) {
        gsl_histogram_increment(h, s);
    }

    gsl_histogram_pdf* pdf = gsl_histogram_pdf_alloc(n_bins);
    gsl_histogram_pdf_init(pdf, h);

    size_t n_mismatches = 0;

    for (size_t i = 1; i < n_bins; ++i) {
        double lower_edge = 0.0;
        double upper_edge = 0.0;

        gsl_histogram_get_range(h, i, &lower_edge, &upper_edge);

        double observed = pdf->sum[i];
        double expected = sampler.cumulative_density(lower_edge);

        if (std::abs(expected - observed) > tolerance) {
            ++n_mismatches;
            WARN("bin " << i << ": "
                 << observed << " != " << expected
                 << " +/- " << tolerance
                 << " (" << expected - observed << ")");
        }
    }

    if (n_mismatches > 0) {
        WARN("* " << n_mismatches << "/" << n_bins
             << " cumulative densities exceed tolerance");
    }

    gsl_histogram_pdf_free(pdf);
    gsl_histogram_free(h);

    return n_mismatches == 0;
}

TEST_CASE("test_samples", "Test sample distribution")
{
    const size_t n_samples = 1000000;
    const size_t n_bins = 1000;
    const double tolerance = 1e-3;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    double lambda = 0.1;

    SECTION("in regime 1") {
        lcfit::rejection_sampler sampler(rng, REGIME_1, lambda);

        std::vector<double> samples = sampler.sample_n(n_samples);
        CHECK(test_sample_distribution(samples, sampler, n_bins, tolerance));
    }

    SECTION("in regime 2") {
        lcfit::rejection_sampler sampler(rng, REGIME_2, lambda);

        std::vector<double> samples = sampler.sample_n(n_samples);
        CHECK(test_sample_distribution(samples, sampler, n_bins, tolerance));
    }

    SECTION("in regime 3") {
        lcfit::rejection_sampler sampler(rng, REGIME_3, lambda);

        std::vector<double> samples = sampler.sample_n(n_samples);
        CHECK(test_sample_distribution(samples, sampler, n_bins, tolerance));
    }

    SECTION("in regime 4") {
        lcfit::rejection_sampler sampler(rng, REGIME_4, lambda);

        std::vector<double> samples = sampler.sample_n(n_samples);
        CHECK(test_sample_distribution(samples, sampler, n_bins, tolerance));
    }

    gsl_rng_free(rng);
}
