#include "catch.hpp"

#include <iostream>
#include "lcfit.h"
#include "lcfit2.h"

double lcfit4_d2f_t(const double t, const bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double r = model->r;
    const double b = model->b;

    const double theta = exp(r * (t + b));

    const double d2 = ((c - m)*pow(r, 2)*pow(theta, 3) - 2*(c + m)*pow(r, 2)*pow(theta, 2) + (c - m)*pow(r, 2)*theta)/(pow(theta, 4) - 2*pow(theta, 2) + 1);

    return d2;
}

TEST_CASE("simple normalized curves are fitted correctly", "[simple_curves]") {
    SECTION("with points sampled from a parabola") {
        const double t0 = 0.1;

        auto f = [t0](double t) {
            return -pow(t - t0, 2.0) - 1000.0;
        };

        std::vector<double> t{0.0, t0, 2*t0};
        std::vector<double> norm_lnl(t.size());

        for (size_t i = 0; i < t.size(); ++i) {
            norm_lnl[i] = f(t[i]) - f(t0);
        }

        const double d1 = 0.0;
        const double d2 = -2.0;

        lcfit2_bsm_t fit_model = {1100, 800, t0, d1, d2};

        lcfit2n_fit(t.size(), t.data(), norm_lnl.data(), &fit_model);

        for (size_t i = 0; i < t.size(); ++i) {
            // epsilon relaxed since the parabola does not have a
            // horizontal asymptote; lcfit2 will not fit it exactly
            REQUIRE(lcfit2_norm_lnl(t[i], &fit_model) == Approx(norm_lnl[i]).epsilon(1e-3));
        }
    }

    SECTION("with points sampled from an lcfit4 log-likelihood function") {
        bsm_t true_model = {1200.0, 800.0, 2.0, 0.5};

        auto f = [&true_model](double t) {
            return lcfit_bsm_log_like(t, &true_model);
        };

        const double t0 = lcfit_bsm_ml_t(&true_model);

        std::vector<double> t{0.0, t0, lcfit_bsm_infl_t(&true_model), 10.0};
        std::vector<double> norm_lnl(t.size());

        for (size_t i = 0; i < t.size(); ++i) {
            norm_lnl[i] = f(t[i]) - f(t0);
        }

        const double d1 = 0.0;
        const double d2 = lcfit4_d2f_t(t0, &true_model);

        lcfit2_bsm_t fit_model = {1100.0, 800.0, t0, d1, d2};

        lcfit2n_fit(t.size(), t.data(), norm_lnl.data(), &fit_model);

        for (size_t i = 0; i < t.size(); ++i) {
            REQUIRE(lcfit2_norm_lnl(t[i], &fit_model) == Approx(norm_lnl[i]));
        }

        bsm_t fit_model4;
        lcfit2_to_lcfit4(&fit_model, &fit_model4);

        REQUIRE(fit_model4.c == Approx(true_model.c));
        REQUIRE(fit_model4.m == Approx(true_model.m));
        REQUIRE(fit_model4.r == Approx(true_model.r));
        REQUIRE(fit_model4.b == Approx(true_model.b));
    }
}

TEST_CASE("automatic fitting of a function works correctly", "[lcfit2_fit_auto]") {
    SECTION("with an lcfit4 log-likelihood function") {
        bsm_t true_model = {1200.0, 800.0, 2.0, 0.5};

        const double t0 = lcfit_bsm_ml_t(&true_model);
        const double d1 = 0.0;
        const double d2 = lcfit4_d2f_t(t0, &true_model);

        lcfit2_bsm_t fit_model = {1100.0, 800.0, t0, d1, d2};

        const double min_t = 0.0;
        const double max_t = 10.0;
        const double alpha = 0.0;

        double (*f)(double, void*) = reinterpret_cast<double (*)(double, void*)>(&lcfit_bsm_log_like);
        lcfit2_fit_auto(f, &true_model, &fit_model, min_t, max_t, alpha);

        bsm_t fit_model4;
        lcfit2_to_lcfit4(&fit_model, &fit_model4);

        REQUIRE(fit_model4.c == Approx(true_model.c));
        REQUIRE(fit_model4.m == Approx(true_model.m));
        REQUIRE(fit_model4.r == Approx(true_model.r));
        REQUIRE(fit_model4.b == Approx(true_model.b));
    }
}
