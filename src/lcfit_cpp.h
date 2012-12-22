#ifndef LCFIT_CPP_H
#define LCFIT_CPP_H

#include "lcfit.h"

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

struct LCFitResult
{
    std::vector<Point> evaluated_points;
    bsm_t model_fit;
};

Monotonicity monotonicity(const std::vector<Point>&);
std::vector<Point> select_points(std::function<double(double)>, const std::vector<double>& v, const size_t max_points=8);
std::vector<Point> retain_top(const std::vector<Point>&, const size_t);
LCFitResult fit_bsm_log_likelihood(std::function<double(double)>, const bsm_t&, const std::vector<double>&, const size_t max_points=8);

} // lcfit

#endif // LCFIT_CPP_H
