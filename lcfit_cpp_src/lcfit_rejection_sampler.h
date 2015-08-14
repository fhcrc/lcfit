/**
 * \file lcfit_rejection_sampler.h
 * \brief lcfit C++ rejection sampler
 *
 * This file provides a rejection sampler for sampling under the likelihood
 * curve estimated by lcfit.
 */

#ifndef LCFIT_REJECTION_SAMPLER_H
#define LCFIT_REJECTION_SAMPLER_H

#include <utility>
#include <gsl/gsl_rng.h>
#include "lcfit_cpp.h"

namespace lcfit {

class rejection_sampler {
private:
    gsl_rng* rng_;
    lcfit::LCFitResult fit_result_;

    double ml_t_;
    double ml_ll_;

    double t_min_;
    double t_max_;

    double log_auc_;

    mutable size_t n_trials_;
    mutable size_t n_accepts_;

public:
    rejection_sampler(gsl_rng* rng, const lcfit::LCFitResult& fitResult);
    virtual ~rejection_sampler();

    double sample() const;

    double relative_log_likelihood(double t) const;
    double relative_likelihood(double t) const;

    double log_density(double t) const;
    double density(double t) const;

private:
    const std::pair<double, double> find_bounds_easy() const;
    const std::pair<double, double> find_bounds() const;

    double integrate() const;

    void print_acceptance_rate() const;
};

} // namespace lcfit

#endif // LCFIT_REJECTION_SAMPLER_H
