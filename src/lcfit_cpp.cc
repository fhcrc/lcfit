#include "lcfit_cpp.h"

#include <algorithm>
#include <cassert>
#include <iterator>

using namespace std;

namespace lcfit
{

struct PointsByYRev
{
    inline bool operator()(const Point& p1, const Point& p2)
    {
        return p1.y > p2.y;
    };
};

vector<Point> select_points(std::function<double(double)> log_like, const vector<double>& starting_pts, const size_t max_points=8)
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

vector<Point> retain_top(const vector<Point>& points, const size_t n)
{
    if(n >= points.size()) return points;

    vector<Point> sorted = points;
    // Place the
    partial_sort(begin(sorted), begin(sorted) + n, end(sorted), PointsByYRev());
    sorted.resize(n);
    return sorted;
}


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
    assert(false);
}

}

