#include "lcfit_cpp.h"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <iostream>
#include <stdexcept>

using namespace std;

namespace lcfit
{

struct point_by_y_desc
{
    inline bool operator()(const Point& p1, const Point& p2)
    {
        return p1.y > p2.y;
    };
};

struct point_by_x
{
    inline bool operator()(const Point& p1, const Point& p2)
    {
        return p1.x < p2.x;
    };
};

/// Select points for use with lcfit
/// \param log_like Function returning the actual log-likelihood of a branch length
/// \param starting_pts Initial points to sample. {0.1,0.15,0.5} has given good results. More points will be added ot
/// ensure that the function is non-monotonic on the interval.
/// \param max_points Maximum number of points to sample. Passed to lcfit::select_points.
vector<Point> select_points(std::function<double(double)> log_like, const std::vector<double>& starting_pts, const size_t max_points)
{
        vector<Point> points;

        // Try starting points
        for(const double & d : starting_pts) {
            points.push_back({d, log_like(d)});
        }

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

            assert(is_sorted(points.begin(), points.end(),
                   [](const Point& p1, const Point& p2) -> bool { return p1.x < p2.x; }));
        } while(points.size() <= max_points && c != Monotonicity::NON_MONOTONIC);

        return points;
}

/// Select the top \c n points from \c points by log-likelihood
/// \returns The top \c n points, ordered by increasing \c x value.
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

/// Classify \c points by monotonicity.
/// \param points Input points, <b>sorted by increasing x-value</b>
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

/// Fit the BSM
/// \param log_like Function returning the actual log-likelihood of a branch length
/// \param init_model Initial model. <i>will be scaled to speed fit</i>
/// \param sample_points Initial points to sample. {0.1,0.15,0.5} has given good results. More points will be added ot
/// ensure that the function is non-monotonic on the interval.
/// \param max_points Maximum number of points to sample. Passed to lcfit::select_points.
LCFitResult fit_bsm_log_likelihood(std::function<double(double)> log_like, const bsm_t& init_model, const std::vector<double>& sample_points, const size_t max_points)
{
    bsm_t model = init_model;

    const vector<Point> points = lcfit::select_points(log_like, sample_points, max_points);
    const Point p = *std::max_element(begin(points), end(points),
            [](const Point& p1, const Point& p2) -> bool { return p1.y > p2.y; });
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

    const int status = lcfit_fit_bsm(t.size(), t.data(), l.data(), &model);
    if(status) throw runtime_error("fit_ll returned: " + std::to_string(status));
    return {points, std::move(model)};
}

}

