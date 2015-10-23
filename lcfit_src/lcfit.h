/**
 * \file lcfit.h
 * \brief lcfit C API
 *
 * This file provides functions for fitting the binary symmetric model
 * to empirical likelihood data using non-linear least squares
 * methods.
 *
 */
#ifndef LCFIT_H
#define LCFIT_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Binary symmetric model parameters. */
typedef struct {
    /** Number of constant sites. */
    double c;
    /** Number of mutated sites. */
    double m;
    /** Rate of mutation. */
    double r;
    /** Branch length offset */
    double b;
} bsm_t;

/** Default initial conditions. */
extern const bsm_t DEFAULT_INIT;

/** Status codes returned by #lcfit_fit_bsm and #lcfit_fit_bsm_weight. */
typedef enum {
    /** Success. */
    LCFIT_SUCCESS = 0,
    /** Exceeded maximum number of iterations without converging. */
    LCFIT_MAXITER = 1,
    /** A non-specific error occurred. */
    LCFIT_ERROR = 2,
    /** Iterations are not making progress toward a solution. */
    LCFIT_ENOPROG = 27,
    /** Cannot reach the tolerance specified for the objective function. */
    LCFIT_ETOLF = 29,
    /** Cannot reach the tolerance specified for the objective function's gradient. */
    LCFIT_ETOLG = 31
} lcfit_status;

/** Codes returned by #lcfit_bsm_regime indicating the model's parameter regime. */
typedef enum {
    /** Regime unknown. */
    LCFIT_REGIME_UNKNOWN,
    /** Regime 1. */
    LCFIT_REGIME_1,
    /** Regime 2. */
    LCFIT_REGIME_2,
    /** Regime 3. */
    LCFIT_REGIME_3,
    /** Regime 4. */
    LCFIT_REGIME_4
} lcfit_regime;

/** Compute the BSM log-likelihood at a given branch length.
 *
 *  In general,
 *
 *  \f[
 *    L(t|c,m,r,b) = c \log\left(\frac{1+e^{-r (t+b)}}{2}\right)+
 *                   m \log\left(\frac{1-e^{-r (t+b)}}{2}\right)
 *  \f]
 *
 *  This function handles two special cases:
 *
 *   -# If the model is in regime 1 and t is 0.0, the function
 *      properly returns \c -INFINITY.
 *
 *   -# If t is \c INFINITY (see #lcfit_bsm_ml_t), the function
 *      returns the limit
 *
 *      \f[
 *        \lim_{t \to \infty} L(t|c,m,r,b) = (c + m) \log(0.5)
 *      \f]
 *
 *      This result enables rejection sampling even when the model is
 *      in regime 4 and the likelihood curve is monotonically
 *      increasing; see lcfit::rejection_sampler::sample for details.
 *
 *  \param[in] t  Branch length.
 *  \param[in] m  Model parameters.
 *
 *  \return The log-likelihood under \c m.
 */
double lcfit_bsm_log_like(double t, const bsm_t* m);

/** Compute the maximum-likelihood branch length for a given model.
 *
 * In general,
 *
 * \f[
 *   \hat{t} = -b + \frac{1}{r} \log \left( \frac{c + m}{c - m} \right)
 * \f]
 *
 * The branch length is constrained to be non-negative; if the above
 * calculation yields a negative value, the function returns 0.0
 * instead.
 *
 * If the model is in regime 4 (i.e., <c>m.c < m.m</c>), where the
 * likelihood curve is monotonically increasing, the function properly
 * returns \c INFINITY for use with #lcfit_bsm_log_like.
 *
 * \param[in] m  Model parameters.
 *
 * \return The maximum-likelihood branch length under \c m.
 */
double lcfit_bsm_ml_t(const bsm_t* m);

/** Compute the positive inflection point for a model in regime 1 or 2.
 *
 * The second derivative of the surrogate function is zero when
 *
 * \f[
 *   t = -b + \frac{1}{r} \log \left( \frac{(\sqrt{c} \pm \sqrt{m})^2}{c - m} \right)
 * \f]
 *
 * This equation has a positive real solution only when the model is
 * in regime 1 or 2. When it exists, the point is an inflection point.
 *
 * \param[in] m  Model parameters.
 *
 * \return The positive inflection point under \c m, or \c NAN if none.
 */
double lcfit_bsm_infl_t(const bsm_t* m);

/** Compute the model parameter gradient at a given branch length.
 *
 * Let \f$u = e^{-r (t + b)}\f$. This function computes the model
 * parameter gradient as the four partial derivatives
 *
 * \f[
 *   \frac{\partial f}{\partial c} = \log \left( \frac{1 + u}{2} \right)
 * \f]
 *
 * \f[
 *   \frac{\partial f}{\partial m} = \log \left( \frac{1 - u}{2} \right)
 * \f]
 *
 * \f[
 *   \frac{\partial f}{\partial r} = (t + b) \left( -c \frac{u}{1 + u} + m \frac{u}{1 - u} \right)
 * \f]
 *
 * \f[
 *   \frac{\partial f}{\partial b} = r \left( -c \frac{u}{1 + u} + m \frac{u}{1 - u} \right)
 * \f]
 *
 *  \param[in]     t     Branch length.
 *  \param[in]     m     Model parameters.
 *  \param[in,out] grad  A pointer to a preallocated four-element array of type
 *                       \c double for storing the model parameter gradient at
 *                       \c t.
 */
void lcfit_bsm_gradient(const double t, const bsm_t* m, double* grad);

/** Determine the parameter regime for a model.
 *
 * \param[in] m  Model parameters.
 *
 * \return An #lcfit_regime code indicating the model regime.
 */
lcfit_regime lcfit_bsm_regime(const bsm_t* m);

/** Compute a scale factor for a given model to intersect with \f$(t, l)\f$.
 *
 * This function computes a scaling parameter for the values \c m.c
 * and \c m.m to obtain log-likelihood value \c l at branch-length \c
 * t, keeping \c m.r and \c m.b fixed.
 *
 * \param[in] t  Branch length.
 * \param[in] l  Log-likelihood.
 * \param[in] m  Model parameters.
 *
 * \return A value \c v, such that the likelihood curve for the model
 *         <c>{c / v, m / v, r, b}</c> intersects with \f$(t, l)\f$
 */
double lcfit_bsm_scale_factor(const double t, const double l, const bsm_t* m);

/**
 * Rescale a given model to intersect with \f$(t, l)\f$.
 *
 * See #lcfit_bsm_scale_factor.
 *
 * \param[in]     t  Branch length.
 * \param[in]     l  Log-likelihood.
 * \param[in,out] m  Model parameters, updated in-place.
 */
void lcfit_bsm_rescale(const double t, const double l, bsm_t* m);

/** Fit a given model to empirical likelihood data with weighting.
 *
 * This function fits the binary symmetric model to weighted empirical
 * likelihood samples using non-linear least squares methods, starting
 * with the initial conditions specified by \c m. Combine
 * #DEFAULT_INIT and #lcfit_bsm_scale_factor for reasonable starting
 * conditions.
 *
 * \param[in]     n  Number of observations in \c t and \c l.
 * \param[in]     t  Branch lengths.
 * \param[in]     l  Log-likelihood value at each \c t.
 * \param[in]     w  Weight for sample point at each \c t.
 * \param[in,out] m  Model parameters, updated in-place.
 *
 * \return An #lcfit_status code, zero for success, non-zero otherwise.
 */
int lcfit_fit_bsm_weight(const size_t n, const double* t, const double* l,
                         const double* w, bsm_t* m, int max_iter);

/** Fit a given model to empirical likelihood data without weighting.
 *
 * This is a convenience function for fitting a model without sample
 * weighting; it simply calls #lcfit_fit_bsm_weight with equal weights
 * for the samples. See that function's documentation for more
 * information.
 *
 * \param[in]     n  Number of observations in \c t and \c l.
 * \param[in]     t  Branch lengths.
 * \param[in]     l  Log-likelihood value at each \c t.
 * \param[in,out] m  Model parameters, updated in-place.
 *
 * \return An #lcfit_status code, zero for success, non-zero otherwise.
  */
int lcfit_fit_bsm(const size_t n, const double* t, const double* l, bsm_t* m,
                  int max_iter);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // LCFIT_H
