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

/// lcfit C++ API
namespace lcfit
{

/// Function monotonicity
enum class Monotonicity
{
    /// Monotonicity unknown
    UNKNOWN = 0,
    /// Monotonically increasing
    MONO_INC,
    /// Monotonically decreasing
    MONO_DEC,
    /// Non-monotonic (points enclose an extremum)
    NON_MONOTONIC
};

/// A sampled point, with x and y values
struct Point
{
    double x, y;

    inline bool operator==(const Point& other) const
    {
        return x == other.x && y == other.y;
    };
};

/// The result of running lcfit
struct LCFitResult
{
    /// The points evaluated to obtain model_fit
    std::vector<Point> evaluated_points;
    /// The fit BSM
    bsm_t model_fit;
};

/// Classify \c points by monotonicity.
///
/// \param points Input points, <b>sorted by increasing x-value</b>
Monotonicity monotonicity(const std::vector<Point>& points);

/// Select points for use with lcfit
///
/// \param log_like Function returning the actual log-likelihood of a branch length
/// \param starting_pts Initial points to sample. {0.1,0.15,0.5} has given good results. More points will be added to
/// ensure that the function is non-monotonic on the interval.
/// \param max_points Maximum number of points to sample. Passed to lcfit::select_points.
std::vector<Point> select_points(std::function<double(double)>, const std::vector<double>& v, const size_t max_points=8);

std::vector<Point> select_points(std::function<double(double)>, const std::vector<Point>&, const size_t max_points=8);

/// Select the top \c n points from \c points by y-value
///
/// \param points Input points
/// \param n Number of points to retain
/// \returns The top \c n points, ordered by increasing \c x value.
std::vector<Point> retain_top(const std::vector<Point>& points, const size_t n);

/// Fit the BSM.
///
/// Given some initial points to sample, generated as <c>(x, log_like(x))</c>:
/// * chooses additional points so that the collection of points is non-monotonic
/// * Runs #lcfit_fit_bsm on the selected points
///
/// \param log_like Function returning the actual log-likelihood of a branch length
/// \param init_model Initial model. <i>will be scaled to speed fit</i>
/// \param sample_points Initial points to sample. <c>{0.1,0.15,0.5}</c> has given good results. More points will be added ot
/// ensure that the function is non-monotonic on the interval.
/// \param max_points Maximum number of points to sample. Passed to lcfit::select_points.
/// \param max_iter Maximum number of iterations to run. Passed to lcfit::fit_bsm_log_likelihood.
/// \returns The fit BSM, along with points evaluated during fitting
LCFitResult fit_bsm_log_likelihood(std::function<double(double)>, const bsm_t&, const std::vector<double>&, const size_t max_points=8, const int max_iter=250);

} // lcfit

#endif // LCFIT_CPP_H
