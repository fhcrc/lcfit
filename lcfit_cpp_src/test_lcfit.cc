#include "catch.hpp"

#include <vector>
#include "lcfit.h"
#include "lcfit_select.h"

TEST_CASE("model values are computed correctly", "[bsm_values]") {
    SECTION("in regime 1") {
        bsm_t m = {10.0, 1.0, 1.0, 0.0};
        REQUIRE(lcfit_bsm_ml_t(&m) == Approx(0.2006707));
        REQUIRE(lcfit_bsm_log_like(0.0, &m) == -INFINITY);
        REQUIRE(lcfit_bsm_log_like(0.2, &m) == Approx(-3.351002));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &m) == Approx(-7.624619));
    }

    SECTION("in regime 2") {
        bsm_t m = {10.0, 1.0, 1.0, 0.1};
        REQUIRE(lcfit_bsm_ml_t(&m) == Approx(0.1006707));
        REQUIRE(lcfit_bsm_log_like(0.0, &m) == Approx(-3.532821));
        REQUIRE(lcfit_bsm_log_like(0.1, &m) == Approx(-3.351002));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &m) == Approx(-7.624619));
    }

    SECTION("in regime 3") {
        bsm_t m = {10.0, 1.0, 1.0, 1.0};
        REQUIRE(lcfit_bsm_ml_t(&m) == 0.0);
        REQUIRE(lcfit_bsm_log_like(0.0, &m) == Approx(-4.950677));
        REQUIRE(lcfit_bsm_log_like(0.1, &m) == Approx(-5.156038));
        REQUIRE(lcfit_bsm_log_like(INFINITY, &m) == Approx(-7.624619));
    }

    SECTION("in regime 4") {
        bsm_t m = {1.0, 10.0, 1.0, 0.1};
        REQUIRE(lcfit_bsm_ml_t(&m) == INFINITY);
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
    log_like_function_t log_like =
            {(double (*)(double, void*)) &lcfit_bsm_log_like, &model};

    SECTION("when the maximum is to the left of the starting points") {
        std::vector<double> t = {0.5, 1.0, 1.1};
        std::vector<point_t> starting_pts = {{t[0], lcfit_bsm_log_like(t[0], &model)},
                                             {t[1], lcfit_bsm_log_like(t[1], &model)},
                                             {t[2], lcfit_bsm_log_like(t[2], &model)}};
        size_t n_pts = starting_pts.size();
        const size_t max_pts = 8;

        point_t* selected_pts = select_points(&log_like, starting_pts.data(),
                                              &n_pts, max_pts);

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
                                              &n_pts, max_pts);

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
                                              &n_pts, max_pts);

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

TEST_CASE("trying to fit a model returns an error",
          "[lcfit_fit_bsm_errors]") {

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

TEST_CASE("estimated maximum likelihood branch length is within tolerance",
          "[estimate_ml_t]") {
    bsm_t true_model = {1200.0, 300.0, 1.0, 0.2}; // ml_t = 0.310826
    const double true_ml_t = lcfit_bsm_ml_t(&true_model);
    log_like_function_t log_like =
            {(double (*)(double, void*)) &lcfit_bsm_log_like, &true_model};

    bsm_t model = {1800.0, 400.0, 1.0, 0.5};
    bool success = false;

    SECTION("when the tolerance is 1e-3") {
        const double tolerance = 1e-3;

        SECTION("and the maximum is enclosed by the starting points") {
            std::vector<double> t = {0.1, 0.5, 1.0, 1.5};

            double result = estimate_ml_t(&log_like, t.data(), t.size(),
                                          tolerance, &model, &success);

            REQUIRE(success == true);
            REQUIRE(result == Approx(true_ml_t).epsilon(tolerance));
        }

        SECTION("and the maximum is to the left of the starting points") {
            std::vector<double> t = {1.0, 1.1, 1.4, 1.5};

            double result = estimate_ml_t(&log_like, t.data(), t.size(),
                                          tolerance, &model, &success);

            REQUIRE(success == true);
            REQUIRE(result == Approx(true_ml_t).epsilon(tolerance));
        }
    }
}
