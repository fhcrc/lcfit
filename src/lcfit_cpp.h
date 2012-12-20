#ifndef LCFIT_CPP_H
#define LCFIT_CPP_H

#include <cstdlib>
#include <functional>
#include <vector>

namespace lcfit
{

enum class Monotonicity
{
    UNKNOWN = 0,
    MONO_INC,
    MONO_DEC,
    NON_MONOTONIC
};

struct Point
{
    double x, y;

    inline bool operator==(const Point& other) const
    {
        return x == other.x && y == other.y;
    };
};

Monotonicity monotonicity(const std::vector<Point>&);
std::vector<Point> select_points(std::function<double(double)>, const std::vector<double>&, const size_t);
std::vector<Point> retain_top(const std::vector<Point>&, const size_t);

} // lcfit

#endif // LCFIT_CPP_H
