#ifndef STS_GSL_H
#define STS_GSL_H

#include <functional>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>

namespace gsl
{
double minimize(const std::function<double(double)> fn,
                double m = 0.5,
                double a = 0,
                double b = 1,
                const int max_iter = 100,
                const double tolerance = 1e-3,
                const gsl_min_fminimizer_type *min_type = gsl_min_fminimizer_brent);

double find_root(const std::function<double(double)> fn,
                 double a = 0,
                 double b = 1,
                 const int max_iter = 100,
                 const double tolerance = 1e-3,
                 const gsl_root_fsolver_type *min_type = gsl_root_fsolver_brent);
} // namespace gsl

#endif // STS_GSL_H
