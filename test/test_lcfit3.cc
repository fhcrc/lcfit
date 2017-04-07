#include "catch.hpp"

#include <iostream>

#include "lcfit.h"
#include "lcfit3.h"

double lcfit4_d1f_t(const double t, const bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double r = model->r;
    const double b = model->b;

    const double theta = exp(r * (t + b));

    const double d1 = -c*r/(theta*(1/theta + 1)) - m*r/(theta*(1/theta - 1));

    return d1;
}

// defined in test_lcfit2.cc
double lcfit4_d2f_t(const double t, const bsm_t* model);

const bsm_t REGIME_3 = {10.0, 1.0, 1.0, 1.0};

TEST_CASE("lcfit3 automatic fitting of a function works correctly", "[lcfit3_fit_auto]") {
    SECTION("with an lcfit4 log-likelihood function") {
        bsm_t true_model = REGIME_3;

        const double c = 1100.0;
        const double m = 800.0;
        const double theta_b = (c + m + 2*sqrt(c*m))/(c - m) + 1.0;

        const double d1 = lcfit4_d1f_t(0.0, &true_model);
        const double d2 = lcfit4_d2f_t(0.0, &true_model);

        lcfit3_bsm_t fit_model = {c, m, theta_b, d1, d2};

        const double min_t = 0.0;
        const double max_t = 10.0;
        const double alpha = 0.0;

        double (*f)(double, void*) = reinterpret_cast<double (*)(double, void*)>(&lcfit_bsm_log_like);
        lcfit3_fit_auto(f, &true_model, &fit_model, min_t, max_t, alpha);

        bsm_t fit_model4;
        lcfit3_to_lcfit4(&fit_model, &fit_model4);

        CHECK(fit_model4.c == Approx(true_model.c));
        CHECK(fit_model4.m == Approx(true_model.m));
        CHECK(fit_model4.r == Approx(true_model.r));
        CHECK(fit_model4.b == Approx(true_model.b));
    }
}
