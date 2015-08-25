/**
 * \file lcfit_rejection_sampler.h
 * \brief lcfit C++ rejection sampler
 *
 * This file provides a rejection sampler for sampling and computing
 * likelihoods and approximate densities under the BSM likelihood
 * curve with an exponential prior.
 */

#ifndef LCFIT_REJECTION_SAMPLER_H
#define LCFIT_REJECTION_SAMPLER_H

#include <vector>
#include <gsl/gsl_rng.h>

#include "lcfit.h"

namespace lcfit {

/**
 * A rejection sampler for sampling and computing likelihoods and
 * approximate densities under a BSM likelihood curve with an
 * exponential prior.
 *
 * Example usage:
 *
 * \code
 * gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
 *
 * bsm_t model = {10.0, 1.0, 1.0, 0.0};
 * double lambda = 0.1;
 *
 * lcfit::rejection_sampler sampler(rng, model, lambda);
 *
 * std::vector<double> samples(1000);
 * std::generate(samples.begin(), samples.end(),
 *               [&sampler]() { return sampler.sample(); });
 * \endcode
 */
class rejection_sampler {
private:
    gsl_rng* rng_;
    bsm_t model_;
    double mu_;

    double ml_t_;
    double ml_ll_;

    mutable double log_auc_;
    mutable bool log_auc_cached_;

  public:
    /**
     * Construct a new sampler given a model and an exponential prior.
     *
     * \param[in,out] rng     GSL random number generator.
     * \param[in]     model   Model parameters.
     * \param[in]     lambda  Rate of exponential prior.
     */
    rejection_sampler(gsl_rng* rng, const bsm_t& model, double lambda);
    virtual ~rejection_sampler() = default;

    /** Generate a sample from the distribution. */
    double sample() const;

    /** Generate multiple samples from the distribution. */
    std::vector<double> sample_n(size_t n) const;

    /** Compute the log-likelihood at a given branch length. */
    double log_likelihood(double t) const;

    /** Compute the likelihood at a given branch length. */
    double likelihood(double t) const;

    /**
     * Compute the approximate log density at a given branch length.
     *
     * The first time any of the `density` functions is called, the
     * normalization constant will be approximated by numerical
     * integration of the likelihood curve and cached for future
     * calls.
     */
    double log_density(double t) const;

    /**
     * Compute the approximate density at a given branch length.
     *
     * The first time any of the `density` functions is called, the
     * normalization constant will be approximated by numerical
     * integration of the likelihood curve and cached for future
     * calls.
     */
    double density(double t) const;

    /** Compute the approximate cumulative density at a given branch length.
     *
     * The first time any of the `density` functions is called, the
     * normalization constant will be approximated by numerical
     * integration of the likelihood curve and cached for future
     * calls.
     */
    double cumulative_density(double t) const;

private:
    /** Compute the approximate integral of the unnormalized posterior. */
    double integrate(double t) const;
};

} // namespace lcfit

#endif // LCFIT_REJECTION_SAMPLER_H
