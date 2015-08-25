/**
 * \file lcfit_rejection_sampler.h
 * \brief lcfit C++ rejection sampler
 *
 * This file provides a rejection sampler for sampling under the likelihood
 * curve estimated by lcfit.
 */

#ifndef LCFIT_REJECTION_SAMPLER_H
#define LCFIT_REJECTION_SAMPLER_H

#include <vector>
#include <gsl/gsl_rng.h>

#include "lcfit.h"

namespace lcfit {

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
    rejection_sampler(gsl_rng* rng, const bsm_t& model, double lambda);
    virtual ~rejection_sampler() = default;

    double sample() const;
    std::vector<double> sample_n(size_t n) const;

    double log_likelihood(double t) const;
    double likelihood(double t) const;

    double log_density(double t) const;
    double density(double t) const;

private:
    double integrate() const;
};

} // namespace lcfit

#endif // LCFIT_REJECTION_SAMPLER_H
