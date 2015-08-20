#include "lcfit_rejection_sampler.h"

#include <cmath>
#include <stdexcept>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "lcfit_cpp.h"

namespace lcfit {

rejection_sampler::rejection_sampler(gsl_rng* rng, const bsm_t& model, double lambda) :
    rng_(rng), model_(model), mu_(1.0 / lambda)
{
    if (!std::isnormal(mu_)) {
        throw std::invalid_argument("invalid exponential rate parameter");
    }

    ml_t_ = lcfit_bsm_ml_t(&model_);

    if (std::isnan(ml_t_)) {
        throw std::runtime_error("lcfit failure: ML branch length is NaN");
    }

    ml_ll_ = lcfit_bsm_log_like(ml_t_, &model_);

    if (!std::isfinite(ml_ll_)) {
        throw std::runtime_error("lcfit failure: non-finite ML log-likelihood");
    }

    log_auc_ = std::log(integrate());

    if (!std::isnormal(log_auc_)) {
        throw std::runtime_error("lcfit failure: invalid integration result");
    }
}

double rejection_sampler::sample() const
{
    double t = 0.0;
    double u = 0.0;
    double f = 0.0;

    do {
        t = gsl_ran_exponential(rng_, mu_);
        u = gsl_rng_uniform(rng_);
        f = std::exp(lcfit_bsm_log_like(t, &model_) - ml_ll_);
    } while (u >= f);

    return t;
}

double rejection_sampler::relative_log_likelihood(double t) const
{
    return lcfit_bsm_log_like(t, &model_)
            + std::log(gsl_ran_exponential_pdf(t, mu_))
            - ml_ll_;
}

double rejection_sampler::relative_likelihood(double t) const
{
    return std::exp(relative_log_likelihood(t));
}

double rejection_sampler::log_density(double t) const
{
    return relative_log_likelihood(t) - log_auc_;
}

double rejection_sampler::density(double t) const
{
    return std::exp(log_density(t));
}

double relative_likelihood_callback(double t, void* data)
{
    return static_cast<const lcfit::rejection_sampler*>(data)->relative_likelihood(t);
}

double rejection_sampler::integrate() const
{
    gsl_function f;
    f.function = &lcfit::relative_likelihood_callback;
    f.params = static_cast<void*>(const_cast<rejection_sampler*>(this));

    double result = 0.0;
    double error = 0.0;

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(100);
    int status = gsl_integration_qagiu(&f, 0.0, 0.0, 1e-5, 100, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);

    if (status) {
        throw std::runtime_error(gsl_strerror(status));
    }

    return result;
}

} // namespace lcfit
