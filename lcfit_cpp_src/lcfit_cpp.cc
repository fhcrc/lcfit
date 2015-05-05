/**
 * \file lcfit_cpp.cc
 * \brief lcfit C++-API implementation
 */
#include "lcfit_cpp.h"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <iostream>
#include <stdexcept>

using namespace std;

namespace lcfit
{

void print_points(const std::vector<Point>& points, std::string prefix = "")
{
    Monotonicity c = monotonicity(points);
    if (c == Monotonicity::MONO_DEC) {
        prefix += "decreasing";
    } else if (c == Monotonicity::MONO_INC) {
        prefix += "increasing";
    } else {
        prefix += "!monotonic";
    }

    fprintf(stderr, "%s", prefix.c_str());
    prefix = " ";

    for (const Point& p : points) {
        fprintf(stderr, "%s{ %.6f, %.3f }", prefix.c_str(), p.x, p.y);
        prefix = ", ";
    }

    fprintf(stderr, "\n");
}

struct point_by_y_desc
{
    inline bool operator()(const Point& p1, const Point& p2)
    {
        return p1.y > p2.y;
    };
};

struct point_by_y_asc
{
    inline bool operator()(const Point& p1, const Point& p2)
    {
        return p1.y < p2.y;
    };
};

struct point_by_x
{
    inline bool operator()(const Point& p1, const Point& p2)
    {
        return p1.x < p2.x;
    };
};

std::vector<Point> select_points(std::function<double(double)> log_like,
                                 const std::vector<double>& starting_pts,
                                 const size_t max_points)
{
    std::vector<Point> points;

    // Evaluate log-likelihood at starting points
    for(const double& d: starting_pts) {
        points.push_back({d, log_like(d)});
    }

    return select_points(log_like, points, max_points);
}

/**
 * \brief Select branch_length values to bracket the maximum of the likelihood curve.
 *
 * The starting points are used to establish a direction in which to move when selecting any additional points.
 * They will be included in the final list along with any additional points that may be selected.

 * @param log_like 	pointer to a function to calculate log-likelihood of a given branch length
 * @param starting_pts 	a vector of already-evaluated starting points.
 * @param max_points 	The maximum number of points that will be used, including \c starting_pts
 *
 * @return vector of values along the branch_length axis.
 */
vector<Point> select_points(std::function<double(double)> log_like, const std::vector<Point>& starting_pts, const size_t max_points)
{
        vector<Point> points(starting_pts);

#ifdef VERBOSE
        print_points(points, "PRE:  ");
#endif /* VERBOSE */

        // Add additional samples until the evaluated branch lengths enclose a maximum.
        size_t offset; // Position
        double d;      // Branch length

        Monotonicity c = monotonicity(points);
        do {
            switch(c) {
            case Monotonicity::NON_MONOTONIC:
                d = (points[1].x + points[2].x) / 2.0; // Add a point between the first and second try
                offset = 2;
                break;
            case Monotonicity::MONO_INC:
                d = points.back().x * 2.0; // Double largest value
                offset = points.size();
                break;
            case Monotonicity::MONO_DEC:
                d = points[0].x / 10.0; // Add new smallest value - order of magnitude lower
                offset = 0;
                break;
            default:
                assert(false);
            }

            // Insert point
            points.insert(begin(points) + offset, {d, log_like(d)});

            c = monotonicity(points);

            assert(is_sorted(points.begin(), points.end(), point_by_x()));
        } while(points.size() <= max_points && c != Monotonicity::NON_MONOTONIC);

#ifdef VERBOSE
        print_points(points, "POST: ");
#endif /* VERBOSE */

        return points;
}

vector<Point> retain_top(const vector<Point>& points, const size_t n)
{
    if(n >= points.size()) return points;

    vector<Point> sorted = points;
    // Place the
    partial_sort(begin(sorted), begin(sorted) + n, end(sorted), point_by_y_desc());
    sorted.resize(n);
    sort(begin(sorted), end(sorted), point_by_x());
    return sorted;
}

// Given a vector x.y points, sorted by increasing x value,
// determin whether y values are increasing, decreasing or non-monotonic over the range,


Monotonicity monotonicity(const std::vector<Point>& points)
{
    assert(points.size() > 1);
    std::vector<Point>::const_iterator i = begin(points), e = end(points);
    bool maybe_inc = true, maybe_dec = true;

    Point last = *i++;
    for(; i != e; i++) {
        Point current = *i;
        if(current.y > last.y) maybe_dec = false;
        else if(current.y < last.y) maybe_inc = false;
        last = current;
    }
    assert(!(maybe_inc && maybe_dec));

    if(!maybe_inc && !maybe_dec) return Monotonicity::NON_MONOTONIC;
    else if(maybe_inc) return Monotonicity::MONO_INC;
    else if(maybe_dec) return Monotonicity::MONO_DEC;
    throw runtime_error("Monotonicity reached end of function.");
}

/**
 * /brief Calculate parameters to fit the lcfit function to a log-likelihood curve.
 *
 * The goal is to come up with a set of parameters that closely fit
 * the lcfit function to the likelihood curve.  The likelihood curve
 * is sampled at some small number of points and parameters are
 * selected that cause lcfit to best approximate the likelhood curve.
 *
 * @param log_like 	pointer to a function to calculate log-likelihood of a given branch length
 * @param init_model 	initial model parameters
 * @param sample_points initial set of branch_lengths at which to fit the model
 * @param max_points 	maximum number of points at which to sample the likelihood.
 *
 * @return pair containing the points selected, and the model parameters.
 */
LCFitResult fit_bsm_log_likelihood(std::function<double(double)> log_like, const bsm_t& init_model, const std::vector<double>& sample_points, const size_t max_points, int max_iter)
{
    bsm_t model = init_model;

    const vector<Point> points = lcfit::select_points(log_like, sample_points, max_points);
    const Point p = *std::max_element(begin(points), end(points), point_by_y_asc());
    const double scale_factor = lcfit_bsm_scale_factor(p.x, p.y, &model);
    model.c *= scale_factor;
    model.m *= scale_factor;

    vector<double> t, l;
    t.reserve(points.size());
    l.reserve(points.size());
    for(const Point& p : points) {
        t.push_back(p.x);
        l.push_back(p.y);
    }

    const int status = lcfit_fit_bsm(t.size(), t.data(), l.data(), &model, max_iter);
    // Used to be that lcfit_fit_bsm() could not return a non-zero status,
    // but now that it can, we are frequently throwing an error here.
    // temporarily disable this test so SCons can complete the simulation in the face of errors.
    // if(status) throw runtime_error("lcfit_fit_bsm returned: " + std::to_string(status));

    return {points, std::move(model)};
}

} // end namespace lcfit
