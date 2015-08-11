#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include "lcfit.h"

static void nrerror(char error_text[])
{
    fprintf(stderr, "nr error: %s\n", error_text);
    _exit(1);
}

/* From "Numerical Recipes in C", 2e, p 214 */
/* Returns the value ln[gamma(xx)] for xx > 0 */
static float gammln(float xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

/* From "Numerical Recipes in C", 2e, p 215 */
/* Returns ln(n!) */
static float factln(int n)
{
    static float a[101];

    if (n < 0) nrerror("Negative factorial in routine factln");
    if (n <= 1) return 0.0;
    if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
    else return gammln(n+1.0);
}

/*
  def comb(n, r):
  return factorial(n) // factorial(r) // factorial(n-r)
*/

/* From "Numerical Recipes in C", 2e, p 215 */
/* Returns the binomial coefficient nCr(n, k) as a floating-point number */
static float bico(int n, int k)
{
    return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

static float bicoln(int n, int k)
{
    return factln(n)-factln(k)-factln(n-k);
}


/* From "Numerical Recipes in C", 2e, p 216 */
/* Returns the value of the beta function B(z, w) */
static float beta(float z, float w)
{
    return exp(gammln(z)+gammln(w)-gammln(z+w));
}

/*
def K(c,m,l,r, x):
  s=0
  for i in range(0,c+1):
    for j in range(0,m+1):
      if ((j%2)==0):
        s=s+ 1.0*comb(c, i)*comb(m, j)/(i+j+l/r)*x**(i+j)
      else:
        s=s- 1.0*comb(c, i)*comb(m, j)/(i+j+l/r)*x**(i+j)

  s=s*x**(1.0*l/r)
  return s
*/

// TODO: lock down types and casts
/* Inversion sampling helper function */
double lcfit_k_exp_prior(const bsm_t* model, double lambda, double x)
{
    if (x == 0.0) {
        return 0.0;
    }

    // round or floor?
    int c = (int)(0.5 + model->c);
    int m = (int)(0.5 + model->m);

    double lr = lambda / model->r;

    double s = 0.0;

    // TODO: vectorize?
    for (int i = 0; i <= c; ++i) {
        for (int j = 0; j <= m; ++j) {
            // TODO: will need to be done in log space or with GMP
            double term = bico(c, i) * bico(m, j) / (i + j + lr) * pow(x, i + j);

            if (j % 2 == 0) {
                s += term;
            } else {
                s -= term;
            }
        }
    }

    s *= pow(x, lr);

    return s;
}

/*
def inv(c,m,l,r,u):
  a=0
  b=1
  x=0.5
  al=K(c,m,l,r,x)
  while (np.fabs(al-u) > 1e-6):
    x=1.0*(a+b)/2
    al= K(c,m,l,r,x)
    if (al>u):
      b=x
    else:
      a=x

  return x
*/

// TODO: lock down types and casts
double lcfit_inv_exp_prior(const bsm_t* model, double lambda, double u)
{
    double a = 0.0;
    double b = 1.0;
    double x = 0.5;
    double al = lcfit_k_exp_prior(model, lambda, x);

    // bisect until K(x) ~= u
    while (fabs(al - u) > 1e-6) {
        if (al > u) {
            b = x;
        } else {
            a = x;
        }

        x = (a + b) / 2.0;
        al = lcfit_k_exp_prior(model, lambda, x);
    }

    return x;
}

/* Usage:

   size_t N = 1000;
   double* a = calloc(N, sizeof(double));

   double k0 = lcfit_k_exp_prior(model, lambda, exp(-1.0 * model->r * model->b));
   for (size_t i = 0; i < N; ++i) {
       double u = unifrnd(0, 1);

       double x = lcfit_inv_exp_prior(model, lambda, (1 - u) * K0);
       double t = -1.0 / model->r * log(x) - b;

       a[i] = t;
   }

*/
