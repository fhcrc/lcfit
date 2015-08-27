#include "lcfit_rejection_sampler.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "lcfit.h"

namespace lcfit {

rejection_sampler::rejection_sampler(gsl_rng* rng, const bsm_t& model, double lambda) :
    rng_(rng), model_(model), mu_(1.0 / lambda), log_auc_cached_(false)
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
}

/**
 * Rejection sampling generates samples from an arbitrary distribution
 * \f$ f(x) \f$ using a proposal distribution \f$ g(x) \f$ subject
 * only to the constraint that \f$ f(x) \leq c g(x) \f$ for some
 * constant \f$ c > 0 \f$. Let \f$ f(t) \f$ be the unnormalized
 * posterior on branch lengths given an exponential prior,
 *
 * \f[
 *   f(t) = \lambda e^{- \lambda t} \ell(t ~|~ \theta)
 * \f]
 *
 * where \f$ \ell(t ~|~ \theta) \f$ is the likelihood function for the
 * binary symmetric model parameterized by \f$ \theta \f$. Let \f$
 * g(t) \f$ be the pdf of the exponential distribution with rate \f$
 * \lambda \f$,
 *
 * \f[
 *   g(t) = \lambda e^{- \lambda t}
 * \f]
 *
 * Clearly the ratio \f$ f(t) / g(t) = \ell(t ~|~ \theta) \f$, so we
 * choose \f$ c \f$ to be the maximum likelihood value
 *
 * \f[
 *   c = \ell(\hat{t} ~|~ \theta)
 * \f]
 *
 * where \f$ \hat{t} \f$ is the BSM maximum-likelihood branch length
 * and can be computed directly. Then we have the ratio
 *
 * \f[
 *   \frac{f(t)}{c g(t)} =
 *     \frac{\ell(t ~|~ \theta)}{\ell(\hat{t} ~|~ \theta)} \leq 1
 * \f]
 *
 * which satisfies the requirement for rejection sampling.
 *
 * The procedure for generating a sample from the distribution begins
 * by drawing a branch length \f$ t \f$ from the exponential
 * distribution with rate \f$ \lambda \f$ and a value \f$ u \f$ from
 * the uniform distribution over \f$ (0, 1] \f$. If
 *
 * \f[
 *   u \leq \frac{f(t)}{c g(t)} = \frac{\ell(t ~|~ \theta)}
 *                                     {\ell(\hat{t} ~|~ \theta)}
 * \f]
 *
 * the sample is accepted; otherwise, the sample is rejected and the
 * procedure is repeated.
 *
 * In addition to simplifying the calculation, eliminating the prior
 * \f$ g(t) \f$ from the acceptance calculation allows sampling from
 * the distribution even when the BSM maximum likelihood branch length
 * is infinite (i.e., regime 4), since the asymptotic maximum
 * likelihood can still be calculated. See \ref lcfit_bsm_log_like for
 * details.
 */
double rejection_sampler::sample() const
{
    double t = 0.0;
    double u = 0.0;
    double f = 0.0;

    do {
        t = gsl_ran_exponential(rng_, mu_);
        u = 1.0 - gsl_rng_uniform(rng_); // 1 - [0, 1) = (0, 1]
        f = std::exp(lcfit_bsm_log_like(t, &model_) - ml_ll_);
    } while (u > f);

    return t;
}

std::vector<double> rejection_sampler::sample_n(size_t n) const
{
    std::vector<double> samples(n);

    std::generate(samples.begin(), samples.end(),
                  [this]() { return sample(); });

    return samples;
}

double rejection_sampler::log_likelihood(double t) const
{
    return lcfit_bsm_log_like(t, &model_)
            + std::log(gsl_ran_exponential_pdf(t, mu_))
            - ml_ll_;
}

double rejection_sampler::likelihood(double t) const
{
    return std::exp(log_likelihood(t));
}

double rejection_sampler::log_density(double t) const
{
    if (!log_auc_cached_) {
        log_auc_ = std::log(integrate(INFINITY));
        log_auc_cached_ = true;
    }

    return log_likelihood(t) - log_auc_;
}

double rejection_sampler::density(double t) const
{
    return std::exp(log_density(t));
}

double rejection_sampler::cumulative_density(double t) const
{
    if (!log_auc_cached_) {
        log_auc_ = std::log(integrate(INFINITY));
        log_auc_cached_ = true;
    }

    return integrate(t) / std::exp(log_auc_);
}

double likelihood_callback(double t, void* data)
{
    return static_cast<const lcfit::rejection_sampler*>(data)->likelihood(t);
}

double rejection_sampler::integrate(double t) const
{
    gsl_function f;
    f.function = &lcfit::likelihood_callback;
    f.params = static_cast<void*>(const_cast<rejection_sampler*>(this));

    double result = 0.0;
    double error = 0.0;

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(100);
    int status;

    if (t == INFINITY) {
        status = gsl_integration_qagiu(&f, 0.0, 0.0, 1e-5, 100,
                                       workspace, &result, &error);
    } else {
        status = gsl_integration_qag(&f, 0.0, t, 0.0, 1e-5, 100, GSL_INTEG_GAUSS21,
                                     workspace, &result, &error);
    }

    gsl_integration_workspace_free(workspace);

    if (status) {
        throw std::runtime_error(gsl_strerror(status));
    }

    if (!std::isnormal(result)) {
        throw std::runtime_error("lcfit failure: invalid integration result");
    }

    return result;
}

} // namespace lcfit
