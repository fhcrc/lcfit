/**
 * \file lcfit_cpp.h
 * \brief lcfit C++ API
 *
 * This file provides more sophisticated methods of choosing points and running lcfit.
 */

#ifndef LCFIT_CPP_H
#define LCFIT_CPP_H

#include "lcfit.h"

#include <cstdlib>
#include <functional>
#include <vector>

/** lcfit C++ API */
namespace lcfit
{

/** Function monotonicity. */
enum class Monotonicity
{
    /** Curve monotonicity is unknown. */
    UNKNOWN = 0,
    /** Points are monotonically increasing. */
    MONO_INC,
    /** Points are monotonically decreasing. */
    MONO_DEC,
    /// Non-monotonic (points enclose an extremum).
    NON_MONOTONIC
};

/** An (x, y) pair. */
struct Point
{
    double x, y;

    inline bool operator==(const Point& other) const
    {
        return x == other.x && y == other.y;
    };
};

/** Results obtained from lcfit. */
struct LCFitResult
{
    /** The points at which the empirical likelihood was evaluated. */
    std::vector<Point> evaluated_points;
    /** The fitted model parameters. */
    bsm_t model_fit;
};

/** Classify a \c std::vector of \ref Point ordered by increasing branch length. */
Monotonicity monotonicity(const std::vector<Point>& points);

/**
 * Select points to use in fitting the BSM for a given log-likelihood function.
 *
 * If the starting points already enclose a maximum, they are returned
 * as is. Otherwise, additional points are added to the left or right
 * of existing points until a maximum is enclosed or the maximum
 * number of points is reached.
 *
 * This function evaluates the log-likelihood at the starting points
 * before calling the overload
 * \ref select_points(std::function<double(double)>, const std::vector<Point>&, const size_t).
 *
 * \param[in]     log_like      Log-likelihood function.
 * \param[in]     starting_pts  Initial branch lengths, in ascending order.
 * \param[in]     max_points    Maximum number of points to try, including those
 *                              in \c starting_pts.
 *
 * \return A \c std::vector of \ref Point enclosing a maximum.
 */
std::vector<Point> select_points(std::function<double(double)>,
                                 const std::vector<double>& v,
                                 const size_t max_points=8);

/**
 * Select points to use in fitting the BSM for a given log-likelihood function.
 *
 * If the starting points already enclose a maximum, they are returned
 * as is. Otherwise, additional points are added to the left or right
 * of existing points until a maximum is enclosed or the maximum
 * number of points is reached.
 *
 * \param[in]     log_like      Log-likelihood function.
 * \param[in]     starting_pts  Already-evaluated sample points, sorted by
 *                              increasing branch length.
 * \param[in]     max_points    Maximum number of points to try, including those
 *                              in \c starting_pts.
 *
 * \return A \c std::vector of \ref Point enclosing a maximum.
 */
std::vector<Point> select_points(std::function<double(double)>,
                                 const std::vector<Point>&,
                                 const size_t max_points=8);

/**
 * Select the top \c n from a \c std::vector of \ref Point by \c y.
 *
 * \param[in] points  A \c std::vector of \ref Point ordered by increasing branch length.
 * \param[in] n       Number of points to return.
 *
 * \return The top \c n points, ordered by increasing \c x.
 */
std::vector<Point> retain_top(const std::vector<Point>& points, const size_t n);

/** Fit a given model to empirical likelihood data with weighting.
 *
 * This function first evaluates the empirical log-likelihood function
 * at the starting branch lengths, and then calls
 * #lcfit::select_points to obtain a non-monotonic interval if the
 * starting points do not enclose an extremum. The function then fits
 * the binary symmetric model to the empirical likelihood samples
 * using non-linear least squares methods, starting with the initial
 * conditions specified by \c init_model.
 *
 * \param[in] log_like       Empirical log-likelihood function.
 * \param[in] init_model     Initial model parameters.
 * \param[in] sample_points  Starting branch lengths.
 * \param[in] max_points     Maximum number of points to sample.
 * \param[in] max_iter       Maximum number of fit iterations.
 *
 * \return The fitted model parameters and points evaluated during fitting.
 */
LCFitResult fit_bsm_log_likelihood(std::function<double(double)>, const bsm_t&,
                                   const std::vector<double>&,
                                   const size_t max_points=8,
                                   const int max_iter=250);

} // namespace lcfit

#endif // LCFIT_CPP_H
