/**
 * \file lcfit_select.h
 * \brief lcfit sample point selection and ML estimation
 *
 * This file provides functions for selecting sample points on an
 * empirical likelihood curve that enclose a maximum and for iterative
 * fitting and maximum-likelihood branch length estimation.
 */
#ifndef LCFIT_SELECT_H
#define LCFIT_SELECT_H
#include <stdlib.h>
#include <stdbool.h>

#include "lcfit.h"

#ifdef __cplusplus
extern "C" {
#endif

/** A (branch length, log-likelihood) pair. */
typedef struct
{
    /** Branch length. */
    double t;
    /** Log-likelihood at \c t. */
    double ll;
} point_t;

/** A log-likelihood function and additional arguments. */
typedef struct
{
    /** Log-likelihood for a given branch length. */
    double (*fn)(double, void*);
    /** Additional arguments to pass to \c fn. */
    void* args;

} log_like_function_t;


/** Classifications returned by #classify_curve. */
typedef enum {
    /** Curve monotonicity is unknown. */
    CRV_UNKNOWN = 0,
    /** Points are monotonically increasing. */
    CRV_MONO_INC = 1,
    /** Points are monotonically decreasing. */
    CRV_MONO_DEC = 2,
    /** Points enclose a minimum. */
    CRV_ENC_MINIMA = 3,
    /** Points enclose a maximum. */
    CRV_ENC_MAXIMA = 4
} curve_type_t;

/** Classify an array of #point_t ordered by increasing branch length. */
curve_type_t
classify_curve(const point_t[], const size_t);

/**
 * Select points to use in fitting the BSM for a given log-likelihood function.
 *
 * If the starting points already enclose a maximum, they are returned
 * as is. Otherwise, additional points are added to the left or right
 * of existing points until a maximum is enclosed or the maximum
 * number of points is reached.
 *
 * \param[in]     log_like      Log-likelihood function.
 * \param[in]     starting_pts  Initial branch lengths, in ascending order.
 * \param[in,out] num_pts       Number of points in \c starting_points;
 *                              modified in-place to contain the number of
 *                              points returned.
 * \param[in]     max_pts       Maximum number of points to try, including those
 *                              in \c starting_pts.
 * \param[in]     min_t         Lower bound on branch length.
 * \param[in]     max_t         Upper bound on branch length.
 *
 * \return An array of #point_t enclosing a maximum, or \c NULL if unable to do so.
 */
point_t*
select_points(log_like_function_t *log_like, const point_t starting_pts[],
              size_t *num_pts, size_t max_pts, const double t_min,
              const double t_max);

/** Sort an array of #point_t by increasing branch length. */
void
sort_by_t(point_t points[], const size_t n);

/** Sort an array of #point_t by decreasing log-likelihood. */
void
sort_by_like(point_t points[], const size_t n);

/**
 * Estimate the ML branch length of a log-likelihood function using lcfit.
 *
 * In essence, this function uses an iteratively-fit model to guide
 * the search for the empirical log-likelihood function's
 * maximum-likelihood branch length.
 *
 * Starting with the initial model given by \c model, this function
 * iteratively refines the model by fitting it to the supplied
 * log-likelihood function with #lcfit_fit_bsm. At each iteration, the
 * procedure evaluates the empirical function at the fit model's ML
 * branch length, computed using #lcfit_bsm_ml_t, and uses the result
 * to update a best guess for the empirical ML branch length. This
 * procedure repeats until the model's computed ML branch length is
 * within the specified tolerance of the estimated empirical ML branch
 * length.
 *
 * This function is intended to offer an alternative to more general
 * function value optimization procedures, such as Brent's method, in
 * situations where evaluating the empirical likelihood function is a
 * significant bottleneck. This method generally requires fewer
 * likelihood evaluations than Brent's method, especially when high
 * precision is not required and the tolerance can be relaxed.
 *
 * If the return value is not \c NAN but \c success is \c false, the
 * routine exceeded the maximum number of iterations before it
 * converged. In this case, the return value, as well as the updated
 * model parameters, represent the best estimate obtained before
 * termination.
 *
 * \param[in]      log_like   Log-likelihood function.
 * \param[in,out]  t          Starting points at which to evaluate
 *                            log-likelihood; modified to store the the top
 *                            \c n_pts points (by likelihood) used for fitting.
 * \param [in]     n_pts      Number of points in \c t.
 * \param [in]     tolerance  Required fit tolerance.
 * \param [in,out] model      Model parameters, updated in-place.
 * \param [out]    success    A \c bool indicating success or failure.
 * \param [in]       min_t    Lower bound on branch length.
 * \param [in]       max_t    Upper bound on branch length.
 *
 * \return The estimated ML branch length, or \c NAN if an error occurs.
  */
double
estimate_ml_t(log_like_function_t *log_like, double t[],
              const size_t n_pts, const double tolerance, bsm_t* model,
              bool* success, const double min_t, const double max_t);

/**
 * Choose the top \c k points by log-likelihood while maintaining monotonicity.
 *
 * The point with the maximum log-likelihood and the points on either
 * side of it are chosen first. Then <c>k - 3</c> additional points
 * with the highest log-likelihood are chosen.
 *
 * Exposed for testing.
 *
 * \param[in,out] p  Array of #point_t ordered by increasing branch length.
 * \param[in]     n  Number of points in \c p.
 * \param[in]     k  Number of points to return.
 */
void
subset_points(point_t p[], const size_t n, const size_t k);

#ifdef LCFIT_DEBUG
void
lcfit_select_initialize(void);
#endif

#ifdef __cplusplus
} // extern "C"
#endif

#endif
