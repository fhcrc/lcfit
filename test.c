#include <stdlib.h>
#include <stdio.h>
#include "lcfit.h"

/*
Should return:
fit values: 5	1	1
t_hat: 0.405465
*/

int main(void)
{
    const size_t n = 3;
    double t[3] = {0.25, 0.5, 1};
    double l[3] = { -2.787877533411486, -2.721250292026327, -3.05124979115564};
    double x[3] = { 100.0, 2.0, 1.0 };

    int status = fit_ll(n, t, l, x);

    double c = x[0];
    double r = x[1];
    double m = x[2];

    int i;
    printf("fit values: ");
    for(i=0;i<3;i++)
      printf("%g\t", x[i]);
    printf("\n");

    double t_hat = ml_t(c, m, r);

    printf("t_hat: %g\n", t_hat);

    return status;
}

