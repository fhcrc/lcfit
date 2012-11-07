#include <stdlib.h>
#include <stdio.h>
#include "lcfit.h"

int main(void)
{
    const size_t n = 3;
    double t[3] = {0.25, 0.5, 1};
    double l[3] = { -2.787877533411486, -2.721250292026327, -3.05124979115564};
    double x_init[3] = { 100.0, 2.0, 1.0 };

    int status = fit(n, t, l, x_init);

    int i;
    for(i=0;i<3;i++)
      printf("%g\t", x_init[i]);

    printf("\n");

    return status;
}

