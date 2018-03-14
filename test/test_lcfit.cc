#include "catch.hpp"

#include <vector>
#include "lcfit.h"
#include "lcfit_priv.h"
#include "lcfit_select.h"
#include <gsl/gsl_version.h>

std::ostream& operator<<(std::ostream& os, const bsm_t& model)
{
    os << "["
       << model.c << ", "
       << model.m << ", "
       << model.r << ", "
       << model.b << "]";

    return os;
}

double lcfit_lnl_callback(double t, void* data)
{
    const bsm_t* model = reinterpret_cast<const bsm_t*>(data);
    return lcfit_bsm_log_like(t, model);
}

const bsm_t REGIME_1 = {10.0, 1.0, 1.0, 0.0};
const bsm_t REGIME_2 = {10.0, 1.0, 1.0, 0.1};
const bsm_t REGIME_3 = {10.0, 1.0, 1.0, 1.0};
const bsm_t REGIME_4 = {1.0, 10.0, 1.0, 0.1};

const double MIN_BL = 1e-6;
const double MAX_BL = 1e4;

TEST_CASE("model values are computed correctly", "[bsm_values]") {
    SECTION("in regime 1") {
        bsm_t m = REGIME_1;
        REQUIRE(lcfit_bsm_regime(&m) == LCFIT_REGIME_1);
        REQUIRE(lcfit_bsm_ml_t(&m) == Approx(0.2006707));
        REQUIRE(lcfit_bsm_infl_t(&m) == Approx(0.6549003));
        REQUIRE(lcfit_bsm_log_like(0.0, &m) == -INFINITY);
        REQUIRE(lcfit_bsm_log_like(0.2, &m) == Approx(-3.351002));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &m) == Approx(-7.624619));
    }

    SECTION("in regime 2") {
        bsm_t m = REGIME_2;
        REQUIRE(lcfit_bsm_regime(&m) == LCFIT_REGIME_2);
        REQUIRE(lcfit_bsm_ml_t(&m) == Approx(0.1006707));
        REQUIRE(lcfit_bsm_infl_t(&m) == Approx(0.5549003));
        REQUIRE(lcfit_bsm_log_like(0.0, &m) == Approx(-3.532821));
        REQUIRE(lcfit_bsm_log_like(0.1, &m) == Approx(-3.351002));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &m) == Approx(-7.624619));
    }

    SECTION("in regime 3") {
        bsm_t m = REGIME_3;
        REQUIRE(lcfit_bsm_regime(&m) == LCFIT_REGIME_3);
        REQUIRE(lcfit_bsm_ml_t(&m) == 0.0);
        REQUIRE(isnan(lcfit_bsm_infl_t(&m)));
        REQUIRE(lcfit_bsm_log_like(0.0, &m) == Approx(-4.950677));
        REQUIRE(lcfit_bsm_log_like(0.1, &m) == Approx(-5.156038));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &m) == Approx(-7.624619));
    }

    SECTION("in regime 4") {
        bsm_t m = REGIME_4;
        REQUIRE(lcfit_bsm_regime(&m) == LCFIT_REGIME_4);
        REQUIRE(lcfit_bsm_ml_t(&m) == INFINITY);
        REQUIRE(isnan(lcfit_bsm_infl_t(&m)));
        REQUIRE(lcfit_bsm_log_like(0.0, &m) == Approx(-30.50191));
        REQUIRE(lcfit_bsm_log_like(0.1, &m) == Approx(-24.1042));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &m) == Approx(-7.624619));
    }
}

TEST_CASE("curves are classified correctly", "[classify_curve]") {
    SECTION("when points are decreasing") {
        std::vector<point_t> pts = {{0.0, 1.0},
                                    {1.0, 0.8},
                                    {2.0, 0.3},
                                    {3.0, 0.0}};
        REQUIRE(classify_curve(pts.data(), pts.size()) == CRV_MONO_DEC);
    }

    SECTION("when points are eventually decreasing") {
        std::vector<point_t> pts = {{0.0, 1.0},
                                    {1.0, 1.0},
                                    {2.0, 1.0},
                                    {3.0, 0.0}};
        REQUIRE(classify_curve(pts.data(), pts.size()) == CRV_MONO_DEC);
    }

    SECTION("when points are increasing") {
        std::vector<point_t> pts = {{0.0, 1.0},
                                    {1.0, 10.0},
                                    {2.0, 100.0},
                                    {3.0, 1000.0}};
        REQUIRE(classify_curve(pts.data(), pts.size()) == CRV_MONO_INC);
    }

    SECTION("when points are eventually increasing") {
        std::vector<point_t> pts = {{0.0, 0.0},
                                    {1.0, 0.0},
                                    {2.0, 0.0},
                                    {3.0, 1.0}};
        REQUIRE(classify_curve(pts.data(), pts.size()) == CRV_MONO_INC);
    }

    SECTION("when points enclose a maximum") {
        std::vector<point_t> pts = {{0.0, 1.0},
                                    {1.0, 10.0},
                                    {2.0, 100.0},
                                    {3.0, 10.0}};
        REQUIRE(classify_curve(pts.data(), pts.size()) == CRV_ENC_MAXIMA);
    }

    SECTION("when points enclose a minimum") {
        std::vector<point_t> pts = {{0.0, -1.0},
                                    {1.0, -10.0},
                                    {2.0, -100.0},
                                    {3.0, -10.0}};
        REQUIRE(classify_curve(pts.data(), pts.size()) == CRV_ENC_MINIMA);
    }

    SECTION("when points barely enclose a minimum") {
        std::vector<point_t> pts = {{1.0000000000000001e-05, -33047.506890500852},
                                    {0.0001, -33047.506890500852},
                                    {0.001, -33047.506890500852},
                                    {0.01, -33047.506890500852},
                                    {0.10000000000000001, -33047.506890500867},
                                    {0.30000000000000004, -33047.50689050091},
                                    {0.5, -33047.506890500918},
                                    {1, -33047.506890500888}};
        REQUIRE(classify_curve(pts.data(), pts.size()) == CRV_ENC_MINIMA);
    }
}

TEST_CASE("test points are selected properly", "[select_points]") {
    bsm_t model = {1200.0, 300.0, 1.0, 0.2}; // ml_t = 0.310826
    log_like_function_t log_like = {lcfit_lnl_callback, &model};

    SECTION("when the maximum is to the left of the starting points") {
        std::vector<double> t = {0.5, 1.0, 1.1};
        std::vector<point_t> starting_pts = {{t[0], lcfit_bsm_log_like(t[0], &model)},
                                             {t[1], lcfit_bsm_log_like(t[1], &model)},
                                             {t[2], lcfit_bsm_log_like(t[2], &model)}};
        size_t n_pts = starting_pts.size();
        const size_t max_pts = 8;

        point_t* selected_pts = select_points(&log_like, starting_pts.data(),
                                              &n_pts, max_pts, 0.0, INFINITY);

        REQUIRE(n_pts == 4);

        REQUIRE(selected_pts[0].t == Approx(t[0] / 10.0));
        REQUIRE(selected_pts[1].t == Approx(t[0]));
        REQUIRE(selected_pts[2].t == Approx(t[1]));
        REQUIRE(selected_pts[3].t == Approx(t[2]));

        REQUIRE(selected_pts[1].ll > selected_pts[0].ll);
        REQUIRE(selected_pts[1].ll > selected_pts[2].ll);
    }

    SECTION("when the maximum is to the right of the starting points") {
        std::vector<double> t = {0.01, 0.1, 0.2};
        std::vector<point_t> starting_pts = {{t[0], lcfit_bsm_log_like(t[0], &model)},
                                             {t[1], lcfit_bsm_log_like(t[1], &model)},
                                             {t[2], lcfit_bsm_log_like(t[2], &model)}};
        size_t n_pts = starting_pts.size();
        const size_t max_pts = 8;

        point_t* selected_pts = select_points(&log_like, starting_pts.data(),
                                              &n_pts, max_pts, 0.0, INFINITY);

        REQUIRE(n_pts == 5);

        REQUIRE(selected_pts[0].t == Approx(t[0]));
        REQUIRE(selected_pts[1].t == Approx(t[1]));
        REQUIRE(selected_pts[2].t == Approx(t[2]));
        REQUIRE(selected_pts[3].t == Approx(2.0 * t[2]));
        REQUIRE(selected_pts[4].t == Approx(2.0 * 2.0 * t[2]));

        REQUIRE(selected_pts[3].ll > selected_pts[2].ll);
        REQUIRE(selected_pts[3].ll > selected_pts[4].ll);
    }

    SECTION("when the maximum is enclosed by the starting points") {
        std::vector<double> t = {0.1, 0.5, 1.0};
        std::vector<point_t> starting_pts = {{t[0], lcfit_bsm_log_like(t[0], &model)},
                                             {t[1], lcfit_bsm_log_like(t[1], &model)},
                                             {t[2], lcfit_bsm_log_like(t[2], &model)}};
        size_t n_pts = starting_pts.size();
        const size_t max_pts = 8;

        point_t* selected_pts = select_points(&log_like, starting_pts.data(),
                                              &n_pts, max_pts, 0.0, INFINITY);

        REQUIRE (n_pts == 3);

        REQUIRE(selected_pts[0].t == Approx(t[0]));
        REQUIRE(selected_pts[1].t == Approx(t[1]));
        REQUIRE(selected_pts[2].t == Approx(t[2]));

        REQUIRE(selected_pts[1].ll > selected_pts[0].ll);
        REQUIRE(selected_pts[1].ll > selected_pts[2].ll);
    }
}

TEST_CASE("test points are sorted properly", "[sort_points]") {
    auto points_equal = [](const point_t lhs, const point_t rhs) {
        return (lhs.t == Approx(rhs.t)) && (lhs.ll == Approx(rhs.ll));
    };

    SECTION("when sorted by increasing branch length") {
        std::vector<point_t> pts = {{1.0, -4.0},
                                    {0.5, -500.0},
                                    {0.75, -0.75},
                                    {2.0, -28.0}};

        auto points_by_t = [](const point_t lhs, const point_t rhs) {
            return lhs.t < rhs.t;
        };
        std::vector<point_t> sorted_pts = pts;
        std::sort(sorted_pts.begin(), sorted_pts.end(), points_by_t);

        sort_by_t(pts.data(), pts.size());

        REQUIRE(std::equal(pts.begin(), pts.end(), sorted_pts.begin(),
                           points_equal));
    }

    SECTION("when sorted by decreasing likelihood") {
        std::vector<point_t> pts = {{1.0, -4.0},
                                    {0.5, -500.0},
                                    {0.75, -0.75},
                                    {2.0, -28.0}};

        auto points_by_ll_desc = [](const point_t lhs, const point_t rhs) {
            return lhs.ll > rhs.ll;
        };
        std::vector<point_t> sorted_pts = pts;
        std::sort(sorted_pts.begin(), sorted_pts.end(), points_by_ll_desc);

        sort_by_like(pts.data(), pts.size());

        REQUIRE(std::equal(pts.begin(), pts.end(), sorted_pts.begin(),
                           points_equal));
    }
}

TEST_CASE("test point subsets are selected properly", "[subset_points]") {
    auto _test = [](const std::vector<point_t>& expected,
                    std::vector<point_t>& actual,
                    const size_t k) {
        auto points_equal = [](const point_t lhs, const point_t rhs) {
            return (lhs.t == Approx(rhs.t)) && (lhs.ll == Approx(rhs.ll));
        };

        REQUIRE(expected.size() >= k);

        subset_points(actual.data(), actual.size(), k);

        REQUIRE(std::equal(actual.begin(), actual.end(),
                           expected.begin(), points_equal));
    };

    auto order_by = [](const std::vector<point_t>& points,
                       const std::vector<size_t>& indexes) {
        std::vector<point_t> result;
        for (size_t i : indexes) {
            result.push_back(points[i]);
        }
        return result;
    };

    SECTION("when n equals k") {
        std::vector<point_t> pts = {{0.1, 0.4},
                                    {0.2, 0.5},
                                    {0.3, 0.6},
                                    {0.4, 0.1}};

        std::vector<point_t> subset = order_by(pts, {0, 1, 2, 3});
        _test(subset, pts, 4);
    }

    SECTION("when n is less than k") {
        SECTION("and the maximum is at index 1") {
            std::vector<point_t> pts = {{0.1, 0.4},
                                        {0.2, 0.5},
                                        {0.3, 0.3},
                                        {0.4, 0.1}};

            std::vector<point_t> subset = order_by(pts, {0, 1, 2, 3});
            _test(subset, pts, 3);
        }
        SECTION("and the maximum is at index 2") {
            std::vector<point_t> pts = {{0.1, 0.4},
                                        {0.2, 0.5},
                                        {0.3, 0.6},
                                        {0.4, 0.1}};

            std::vector<point_t> subset = order_by(pts, {1, 2, 3, 0});
            _test(subset, pts, 3);
        }
        SECTION("and the maximum is at index 3") {
            std::vector<point_t> pts = {{0.1, 0.4},
                                        {0.2, 0.5},
                                        {0.3, 0.6},
                                        {0.4, 0.8},
                                        {0.45, 0.74},
                                        {0.5, 0.6}};

            std::vector<point_t> subset = order_by(pts, {2, 3, 4, 5, 1, 0});
            _test(subset, pts, 4);
        }
    }
}

TEST_CASE("trying to fit a model returns an error", "[lcfit_fit_bsm_errors]") {
    const int MAX_ITER = 250;

    SECTION("when given fewer than four points to fit") {
        std::vector<double> t = {0.1, 0.5, 1.0};
        std::vector<double> l = {-100.0, -50.0, -90.0};

        bsm_t model = DEFAULT_INIT;
        int status = lcfit_fit_bsm(t.size(), t.data(), l.data(), &model,
                                   MAX_ITER);

        REQUIRE(status == LCFIT_ERROR);
    }
}

TEST_CASE("estimate_ml_t converges to the correct model", "[estimate_ml_t]") {
    const std::vector<double> t = {MIN_BL, 0.1, 0.5, MAX_BL};
    const double tolerance = 1e-3;
    bsm_t fit_model = DEFAULT_INIT;
    bool success = false;

    SECTION("in regime 1") {
        bsm_t true_model = REGIME_1;
        log_like_function_t log_like = {lcfit_lnl_callback, &true_model};

        double fit_ml_t = estimate_ml_t(&log_like, t.data(), t.size(),
                                        tolerance, &fit_model, &success,
                                        MIN_BL, MAX_BL);

        REQUIRE(success == true);
        REQUIRE(fit_ml_t == Approx(lcfit_bsm_ml_t(&true_model)));

        REQUIRE(fit_model.c == Approx(true_model.c));
        REQUIRE(fit_model.m == Approx(true_model.m));
        REQUIRE(fit_model.r == Approx(true_model.r));
        REQUIRE(fit_model.b == Approx(true_model.b));
    }

    SECTION("in regime 2") {
        bsm_t true_model = REGIME_2;
        log_like_function_t log_like = {lcfit_lnl_callback, &true_model};

        double fit_ml_t = estimate_ml_t(&log_like, t.data(), t.size(),
                                        tolerance, &fit_model, &success,
                                        MIN_BL, MAX_BL);

        REQUIRE(success == true);
        REQUIRE(fit_ml_t == Approx(lcfit_bsm_ml_t(&true_model)));

        REQUIRE(fit_model.c == Approx(true_model.c));
        REQUIRE(fit_model.m == Approx(true_model.m));
        REQUIRE(fit_model.r == Approx(true_model.r));
        REQUIRE(fit_model.b == Approx(true_model.b));
    }

    // For this test, estimate_ml_t does not typically converge to the
    // correct model for regimes 3 and 4. The sample points it selects
    // will all be very close to one extreme of the allowable branch
    // length range and thus the empirical curve's behavior at the
    // other extreme is not captured well.
}

TEST_CASE("maxima are bracketed and found properly", "[bracket_and_find]") {
    double min_t = MIN_BL;
    double max_t = MAX_BL;

    SECTION("for curves in regime 1") {
        bsm_t true_model = REGIME_1;
        double true_ml_t = lcfit_bsm_ml_t(&true_model);

        bool success = bracket_maximum(lcfit_lnl_callback, &true_model, &min_t, &max_t);
        REQUIRE(success == true);

        REQUIRE(min_t < true_ml_t);
        REQUIRE(max_t > true_ml_t);

        double guess = (min_t + max_t) / 2.0;

        double est_ml_t = find_maximum(lcfit_lnl_callback, &true_model, guess, min_t, max_t);
        REQUIRE(est_ml_t == Approx(true_ml_t));

        double d1;
        double d2;
        estimate_derivatives(lcfit_lnl_callback, &true_model, est_ml_t, &d1, &d2);

        REQUIRE(d1 == Approx(0.0));
        REQUIRE(d2 < -0.1);
    }

    SECTION("for curves in regime 2") {
        bsm_t true_model = REGIME_2;
        double true_ml_t = lcfit_bsm_ml_t(&true_model);

        bool success = bracket_maximum(lcfit_lnl_callback, &true_model, &min_t, &max_t);
        REQUIRE(success == true);

        REQUIRE(min_t < true_ml_t);
        REQUIRE(max_t > true_ml_t);

        double guess = (min_t + max_t) / 2.0;

        double est_ml_t = find_maximum(lcfit_lnl_callback, &true_model, guess, min_t, max_t);
        REQUIRE(est_ml_t == Approx(true_ml_t));

        double d1;
        double d2;
        estimate_derivatives(lcfit_lnl_callback, &true_model, est_ml_t, &d1, &d2);

        REQUIRE(d1 == Approx(0.0));
        REQUIRE(d2 < -0.1);
    }

    SECTION("for curves in regime 3") {
        bsm_t true_model = REGIME_3;

        bool success = bracket_maximum(lcfit_lnl_callback, &true_model, &min_t, &max_t);
        REQUIRE(success == false);

        double guess = (min_t + max_t) / 2.0;
        REQUIRE(guess == Approx(MIN_BL));
    }

    SECTION("for curves in regime 4") {
        bsm_t true_model = REGIME_4;

        bool success = bracket_maximum(lcfit_lnl_callback, &true_model, &min_t, &max_t);
        REQUIRE(success == false);

        // In regime 4, bracket_maximum will return bounds for which
        // the behavior of the function is nearly indistinguishable
        // (to within the precision of a double) from the asymptotic
        // behavior. As such we test for the proper function behavior
        // rather than testing the bounds themselves.

        double guess = (min_t + max_t) / 2.0;
        REQUIRE(lcfit_lnl_callback(guess, &true_model) == Approx(lcfit_lnl_callback(MAX_BL, &true_model)));

        double d1;
        double d2;
        estimate_derivatives(lcfit_lnl_callback, &true_model, guess, &d1, &d2);

        REQUIRE(d1 == Approx(0.0));
        REQUIRE(d2 == Approx(0.0));
    }
}

double log_callback(double x, void*)
{
    return std::log(x);
}

TEST_CASE("derivatives are estimated properly", "[estimate_derivatives]") {
    double d1;
    double d2;

    SECTION("for log(x), x = 0.1") {
        double x = 0.1;
        estimate_derivatives(&log_callback, NULL, x, &d1, &d2);

        REQUIRE(d1 == Approx(1.0 / x));
        REQUIRE(d2 == Approx(-1.0 / (x * x)));
    }

    SECTION("for log(x), x = 1.0") {
        double x = 1.0;
        estimate_derivatives(&log_callback, NULL, x, &d1, &d2);

        REQUIRE(d1 == Approx(1.0 / x));
        REQUIRE(d2 == Approx(-1.0 / (x * x)));
    }

    SECTION("for log(x), x = 10.0") {
        double x = 10.0;
        estimate_derivatives(&log_callback, NULL, x, &d1, &d2);

        REQUIRE(d1 == Approx(1.0 / x));
        REQUIRE(d2 == Approx(-1.0 / (x * x)));
    }
}

TEST_CASE("lcfit_fit_auto converges to a good model", "[lcfit_fit_auto]") {
    // lcfit2, which lcfit_fit_auto will use under the hood for curves
    // in regimes 1 and 2, doesn't converge to the true regime 1 model
    // {10.0, 1.0, 1.0, 0.0} given the (current DEFAULT_INIT) model
    // {1100.0, 800.0, 2.0, 0.5}. It is unknown why this is the case.
    // lcfit2 actually only uses the c and m values from the initial
    // model, and other values of c and m converge properly, so it is
    // at least known that these model parameters are influencing
    // lcfit2 convergence.
    //
    // In the meantime, we run the regime 1 and 2 tests with a
    // different initial model that has a smaller m value.

    SECTION("in regime 1") {
        bsm_t true_model = REGIME_1;
        double true_ml_t = lcfit_bsm_ml_t(&true_model);
        double true_ml_ll = lcfit_bsm_log_like(true_ml_t, &true_model);

        //bsm_t fit_model = DEFAULT_INIT;
        bsm_t fit_model = {1100.0, 100.0, 2.0, 0.5};

		double fit_ml_t = lcfit_fit_auto(lcfit_lnl_callback, NULL, &true_model, &fit_model, MIN_BL, MAX_BL);
        double fit_ml_ll = lcfit_bsm_log_like(fit_ml_t, &fit_model);

        CAPTURE(fit_model);
        REQUIRE(fit_ml_t == Approx(true_ml_t));

        REQUIRE(fit_model.c == Approx(true_model.c));
        REQUIRE(fit_model.m == Approx(true_model.m));
        REQUIRE(fit_model.r == Approx(true_model.r));
        REQUIRE(fit_model.b == Approx(true_model.b));

        const std::vector<double> ts = {MIN_BL, true_ml_t, MAX_BL};
        for (const double& t : ts) {
            double fit_norm_lnl = lcfit_bsm_log_like(t, &fit_model) - fit_ml_ll;
            double true_norm_lnl = lcfit_bsm_log_like(t, &true_model) - true_ml_ll;

            CAPTURE(t);
            CHECK(fit_norm_lnl == Approx(true_norm_lnl).epsilon(1e-2));
        }
    }

    SECTION("in regime 2") {
        bsm_t true_model = REGIME_2;
        double true_ml_t = lcfit_bsm_ml_t(&true_model);
        double true_ml_ll = lcfit_bsm_log_like(true_ml_t, &true_model);

        //bsm_t fit_model = DEFAULT_INIT;
        bsm_t fit_model = {1100.0, 100.0, 2.0, 0.5};
        double fit_ml_t = lcfit_fit_auto(lcfit_lnl_callback, NULL, &true_model, &fit_model, MIN_BL, MAX_BL);
        double fit_ml_ll = lcfit_bsm_log_like(fit_ml_t, &fit_model);

        CAPTURE(fit_model);
        REQUIRE(fit_ml_t == Approx(true_ml_t));

        REQUIRE(fit_model.c == Approx(true_model.c));
        REQUIRE(fit_model.m == Approx(true_model.m));
        REQUIRE(fit_model.r == Approx(true_model.r));
        REQUIRE(fit_model.b == Approx(true_model.b));

        const std::vector<double> ts = {MIN_BL, true_ml_t, MAX_BL};
        for (const double& t : ts) {
            double fit_norm_lnl = lcfit_bsm_log_like(t, &fit_model) - fit_ml_ll;
            double true_norm_lnl = lcfit_bsm_log_like(t, &true_model) - true_ml_ll;

            CAPTURE(t);
            CHECK(fit_norm_lnl == Approx(true_norm_lnl).epsilon(1e-2));
        }
    }

    SECTION("in regime 3") {
        bsm_t true_model = REGIME_3;
        double true_ml_t = lcfit_bsm_ml_t(&true_model);

        bsm_t fit_model = DEFAULT_INIT;
		double fit_ml_t = lcfit_fit_auto(lcfit_lnl_callback, NULL, &true_model, &fit_model, MIN_BL, MAX_BL);

        CAPTURE(fit_model);
        REQUIRE(fit_ml_t == Approx(true_ml_t));

        const std::vector<double> ts = {MIN_BL, 1e-3, 1e-2, 1e-1, 1.0, 10.0};
        for (const double& t : ts) {
            double fit_lnl = lcfit_bsm_log_like(t, &fit_model);
            double true_lnl = lcfit_bsm_log_like(t, &true_model);

            CAPTURE(t);
            CHECK(fit_lnl == Approx(true_lnl).epsilon(1e-2));
        }
    }
}

TEST_CASE("estimated maximum likelihood branch length is within tolerance", "[ml_t_tolerance]") {
    bsm_t true_model = {1200.0, 300.0, 1.0, 0.2}; // ml_t = 0.310826
    const double true_ml_t = lcfit_bsm_ml_t(&true_model);
    log_like_function_t log_like = {lcfit_lnl_callback, &true_model};

    bsm_t model = {1800.0, 400.0, 1.0, 0.5};
    bool success = false;

    SECTION("when the tolerance is 1e-3") {
        const double tolerance = 1e-3;

        SECTION("and the maximum is enclosed by the starting points") {
            std::vector<double> t = {0.1, 0.5, 1.0, 1.5};

            double result = estimate_ml_t(&log_like, t.data(), t.size(),
                                          tolerance, &model, &success,
                                          MIN_BL, MAX_BL);

            REQUIRE(success == true);
            REQUIRE(result == Approx(true_ml_t).epsilon(tolerance));
        }

        SECTION("and the maximum is to the left of the starting points") {
            std::vector<double> t = {1.0, 1.1, 1.4, 1.5};

            SECTION("and the bounds include the maximum") {
                double result = estimate_ml_t(&log_like, t.data(), t.size(),
                                              tolerance, &model, &success,
                                              MIN_BL, MAX_BL);

                REQUIRE(success == true);
                REQUIRE(result == Approx(true_ml_t).epsilon(tolerance));
            }

            SECTION("and the bounds exclude the maximum") {
                const double t_min = 0.5;
                double result = estimate_ml_t(&log_like, t.data(), t.size(),
                                              tolerance, &model, &success,
                                              t_min, MAX_BL);

                REQUIRE(success == true);
                REQUIRE(result == Approx(t_min).epsilon(tolerance));
            }
        }

        SECTION("and the maximum is to the right of the starting points") {
            std::vector<double> t = {1e-4, 1e-3, 1e-2, 0.1};

            SECTION("and the bounds include the maximum") {
                double result = estimate_ml_t(&log_like, t.data(), t.size(),
                                              tolerance, &model, &success,
                                              MIN_BL, MAX_BL);

                REQUIRE(success == true);
                REQUIRE(result == Approx(true_ml_t).epsilon(tolerance));
            }
        }
    }
}

TEST_CASE("model scale factors are computed correctly", "[lcfit_bsm_scale_factor]") {
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

    REQUIRE(l == Approx(calc_ll));
}

/* Run a single fit, require that it converge, and the residuals decrase from
 * initial conditions */
void
fail_unless_fit_improves(const bsm_t* m, const double t[4], const double l[4])
{
    bsm_t fit = *m;

    int result = lcfit_fit_bsm(4, t, l, &fit, 500);
    REQUIRE(result == 0);

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

void lcfit_bsm_log_like(size_t n, const double* t, double* lnl, const bsm_t* m)
{
    for (size_t i = 0; i < n; ++i) {
        lnl[i] = lcfit_bsm_log_like(t[i], m);
    }
}

TEST_CASE("fitting actually improves fit", "[lcfit_fit_bsm]") {
    const double t[4] = {0.1, 0.2, 0.5, 1.0};
    double l[4];

    SECTION("in regime 1") {
        const bsm_t m = DEFAULT_INIT;
        lcfit_bsm_log_like(4, t, l, &REGIME_1);
        fail_unless_fit_improves(&m, t, l);
    }

    SECTION("in regime 2") {
        const bsm_t m = DEFAULT_INIT;
        lcfit_bsm_log_like(4, t, l, &REGIME_2);
        fail_unless_fit_improves(&m, t, l);
    }

    SECTION("in regime 3") {
        const bsm_t m = DEFAULT_INIT;
        lcfit_bsm_log_like(4, t, l, &REGIME_3);
        fail_unless_fit_improves(&m, t, l);
    }

    SECTION("in regime 4") {
        const bsm_t m = DEFAULT_INIT;
        lcfit_bsm_log_like(4, t, l, &REGIME_4);
        fail_unless_fit_improves(&m, t, l);
    }

    SECTION("using original test case") {
        /* Scaled c and m from scale factor test */
        const bsm_t m = {20739.66, 13826.44, 1.0, 0.5};
        l[0] = -23912.9; l[1] = -23861.9; l[2] = -23804.3; l[3] = -23820.9;
        fail_unless_fit_improves(&m, t, l);
    }
}

TEST_CASE("maximum-likelihood branch lengths are computed properly", "[lcfit_bsm_ml_t]") {
    SECTION("in regime 1") {
        REQUIRE(lcfit_bsm_ml_t(&REGIME_1) == Approx(0.2006707));
    }

    SECTION("in regime 2") {
        REQUIRE(lcfit_bsm_ml_t(&REGIME_2) == Approx(0.1006707));
    }

    SECTION("in regime 3") {
        REQUIRE(lcfit_bsm_ml_t(&REGIME_3) == 0.0);
    }

    SECTION("in regime 4") {
        REQUIRE(lcfit_bsm_ml_t(&REGIME_4) == INFINITY);
    }
}

TEST_CASE("log-likelihoods are computed properly", "[lcfit_bsm_log_like]") {
    SECTION("in regime 1") {
        REQUIRE(lcfit_bsm_log_like(0.0, &REGIME_1) == -INFINITY);
        REQUIRE(lcfit_bsm_log_like(0.1, &REGIME_1) == Approx(-3.532821));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &REGIME_1) == Approx(-7.624619));
    }

    SECTION("in regime 2") {
        REQUIRE(lcfit_bsm_log_like(0.0, &REGIME_2) == Approx(-3.532821));
        REQUIRE(lcfit_bsm_log_like(0.1, &REGIME_2) == Approx(-3.351002));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &REGIME_2) == Approx(-7.624619));
    }

    SECTION("in regime 3") {
        REQUIRE(lcfit_bsm_log_like(0.0, &REGIME_3) == Approx(-4.950677));
        REQUIRE(lcfit_bsm_log_like(0.1, &REGIME_3) == Approx(-5.156038));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &REGIME_2) == Approx(-7.624619));
    }

    SECTION("in regime 4") {
        REQUIRE(lcfit_bsm_log_like(0.0, &REGIME_4) == Approx(-30.50191));
        REQUIRE(lcfit_bsm_log_like(0.1, &REGIME_4) == Approx(-24.1042));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &REGIME_2) == Approx(-7.624619));
    }
}
