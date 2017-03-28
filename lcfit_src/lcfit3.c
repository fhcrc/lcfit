#include "lcfit3.h"

#include <math.h>

double lcfit3_var_r(const lcfit3_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double theta_b = model->theta_b;
    const double f_1 = model->d1;

    const double r = (f_1 * (pow(theta_b, 2.0) - 1)) / ((m - c) * theta_b + m + c);

    return r;
}

double lcfit3_var_q(const lcfit3_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double theta_b = model->theta_b;

    const double q = (c - m) * theta_b - c - m;

    return q;
}

double lcfit3_var_theta(const double t, const lcfit3_bsm_t* model)
{
    const double theta_b = model->theta_b;

    const double r = lcfit3_var_r(model);
    const double theta = theta_b * exp(r * t);

    return theta;
}

void lcfit3n_gradient(const double t, const lcfit3_bsm_t* model, double* grad)
{
    const double c = model->c;
    const double m = model->m;
    const double theta_b = model->theta_b;
    const double f_1 = model->d1;

    const double r = lcfit3_var_r(model);
    const double q = lcfit3_var_q(model);
    const double theta = lcfit3_var_theta(t, model);

    //
    // This is the gradient of the normalized log-likelihood function
    // f(t) - f(0).
    //

    grad[0] = c*r*t*(theta_b - 1)/(q*theta*(1 + 1.0/theta)) + m*r*t*(theta_b - 1)/(q*theta*(-1 + 1.0/theta)) + log(1 + 1.0/theta) - log(1 + 1.0/theta_b);
    grad[1] = -c*r*t*(theta_b + 1)/(q*theta*(1 + 1.0/theta)) - m*r*t*(theta_b + 1)/(q*theta*(-1 + 1.0/theta)) + log(1 - 1/theta) - log(1 - 1/theta_b);
    grad[2] = c*((2*f_1*t*theta_b/q + r*t*(c - m)/q)/theta - 1/(theta*theta_b))/(1 + 1.0/theta) + c/(pow(theta_b, 2)*(1 + 1.0/theta_b)) + m*((2*f_1*t*theta_b/q + r*t*(c - m)/q)/theta - 1/(theta*theta_b))/(-1 + 1.0/theta) + m/(pow(theta_b, 2)*(-1 + 1.0/theta_b));
}

double lcfit3_lnl(const double t, const lcfit3_bsm_t* model)
{
    const double c = model->c;
    const double m = model->m;
    const double theta_b = model->theta_b;

    const double theta = lcfit3_var_theta(t, model);

    const double lnl = c * log(1 + 1/theta) + m * log(1 - 1/theta) - (c + m) * log(2);

    return lnl;
}

double lcfit3_norm_lnl(const double t, const lcfit3_bsm_t* model)
{
    return lcfit3_lnl(t, model) - lcfit3_lnl(0.0, model);
}
