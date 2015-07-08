/**
 * \file lcfit_select.c
 * \brief Point selection implementation
 */
#include "lcfit_select.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

//#ifdef LCFIT_DEBUG
#include <stdio.h>
//#endif

const static size_t MAX_ITERS = 30;

#ifdef LCFIT_DEBUG
static int is_initialized = false;
#endif
static size_t ml_likelihood_calls = 0;
static size_t bracket_likelihood_calls = 0;

#ifdef LCFIT_DEBUG
void show_likelihood_calls(void)
{
    fprintf(stdout, "LCFIT ML Estimation LL calls: %zu\n", ml_likelihood_calls);
    fprintf(stdout, "LCFIT Bracketing LL calls: %zu\n", bracket_likelihood_calls);
}

void lcfit_select_initialize(void)
{
  if(!is_initialized) {
    is_initialized = true;
    atexit(show_likelihood_calls);
  }
}
#endif

/** Default maximum number of points to evaluate in select_points */
static const size_t DEFAULT_MAX_POINTS = 8;

static void
point_ll_minmax(const point_t points[], const size_t n,
                size_t* mini, size_t* maxi)
{
    assert(n > 0 && "No min/max of 0 points");
    *mini = 0;
    *maxi = 0;
    size_t i;
    const point_t* p = points;
    for(i = 1, ++p; i < n; ++i, ++p) {
        if(p->ll > points[*maxi].ll) {
            *maxi = i;
        }
        if(p->ll < points[*mini].ll) {
            *mini = i;
        }
    }
}


curve_type_t
classify_curve(const point_t points[], const size_t n)
{
    size_t mini = n, maxi = n;

#ifndef NDEBUG
    size_t i;
    const point_t *p = points;
    const point_t *last = p++;
    for(i = 1; i < n; ++i, ++p) {
        assert(p->t >= last->t && "Points not sorted!");
        last = p;
    }
#endif

    point_ll_minmax(points, n, &mini, &maxi);
    assert(mini < n);
    assert(maxi < n);

    const size_t end = n - 1;

    if(mini == 0 && maxi == end) {
        return CRV_MONO_INC;
    } else if(mini == end && maxi == 0) {
        return CRV_MONO_DEC;
    } else if(mini != 0 && mini != end && (maxi == 0 || maxi == end)) {
        return CRV_ENC_MINIMA;
    } else if (maxi != 0 && maxi != end && (mini == 0 || mini == end)) {
        return CRV_ENC_MAXIMA;
    }
    return CRV_UNKNOWN;
}

point_t*
select_points(log_like_function_t *log_like, const point_t starting_pts[],
              size_t *num_pts, const size_t max_pts)
{
    point_t* points = malloc(sizeof(point_t) * (max_pts));
    size_t i, n = *num_pts;
    assert(n > 0 && "num_pts must be > 0");

    /* Copy already-evaluated points */
    memcpy(points, starting_pts, sizeof(point_t) * n);
    /* Initialize the rest - for debugging */
    for(i = n; i < max_pts; ++i) {
        points[i].t = -1; points[i].ll = -1;
    }

    /* Add additional samples until the evaluated branch lengths enclose a
     * maximum. */
    size_t offset = 0;  /* Position to maintain sort order */
    double d = 0.0;     /* Branch length */

    curve_type_t c = classify_curve(points, n);
    while (n < max_pts && c != CRV_ENC_MAXIMA) {
        switch(c) {
            case CRV_MONO_INC:
                /* Double largest value */
                d = points[n - 1].t * 2.0;
                offset = n;
                break;
            case CRV_MONO_DEC:
                /* Add new smallest value - order of magnitude lower */
                d = points[0].t / 10.0;
                offset = 0;
                /* shift */
                memmove(points + 1, points, sizeof(point_t) * n);
                break;
            default:
                free(points);
                return NULL;
        }

        const double l = log_like->fn(d, log_like->args);
        points[offset].t = d;
        points[offset].ll = l;
        bracket_likelihood_calls++;

        c = classify_curve(points, ++n);
    }

    *num_pts = n;

    if(n < max_pts)
        return realloc(points, sizeof(point_t) * n);
    else
        return points;
}

static int
dec_like_cmp(const void *p1, const void *p2)
{
    const point_t *pt1 = (const point_t*) p1,
                  *pt2 = (const point_t*) p2;
    if(pt1->ll < pt2->ll) return 1;
    if(pt1->ll == pt2->ll) return 0;
    return -1;
}

static int
t_cmp(const void *p1, const void *p2)
{
    const point_t *pt1 = (const point_t*) p1,
                  *pt2 = (const point_t*) p2;
    if(pt1->t < pt2->t) return -1;
    if(pt1->t == pt2->t) return 0;
    return 1;
}

void
sort_by_like(point_t points[], const size_t n)
{
    qsort(points, n, sizeof(point_t), dec_like_cmp);
}

void
sort_by_t(point_t points[], const size_t n)
{
    qsort(points, n, sizeof(point_t), t_cmp);
}

/*****************/
/* ML estimation */
/*****************/

static inline const point_t*
max_point(const point_t p[], const size_t n)
{
    assert(n > 0);
    const point_t* m = p++;
    size_t i;

    for(i = 1; i < n; ++i) {
        if(p->ll > m->ll) {
            m = p;
        }
        ++p;
    }
    return m;
}

static inline size_t
point_max_index(const point_t p[], const size_t n)
{
    assert(n > 0 && "Cannot take max_index of 0 points");
    if(n == 1) return 0;
    size_t i, max_i = 0;
    for(i = 1; i < n; ++i) {
        if(p[i].ll > p[max_i].ll)
            max_i = i;
    }
    return max_i;
}

void
subset_points(point_t p[], const size_t n, const size_t k)
{
    assert(classify_curve(p, n) == CRV_ENC_MAXIMA);
    if(k == n) return;
    assert(k <= n);
    size_t max_idx = point_max_index(p, n);

    if(max_idx != 1) {
        /* Always keep the points before and after max_idx */
        const size_t n_before = max_idx - 1;
        const size_t s = n_before * sizeof(point_t);
        point_t *buf = malloc(s);
        memcpy(buf, p, s);
        memmove(p, p + max_idx - 1, sizeof(point_t) * 3);
        memcpy(p + 3, buf, s);
        free(buf);
    }
    sort_by_like(p + 3, n - 3);
    sort_by_t(p, k);
}

/* Fill points by evaluating log_like for each value in ts */
static inline void
evaluate_ll(log_like_function_t *log_like, const double *ts,
            const size_t n_pts, point_t *points)
{
    point_t *p = points;
    size_t i;
    for(i = 0; i < n_pts; ++i, ++p) {
        p->t = *ts++;
        p->ll = log_like->fn(p->t, log_like->args);
        ml_likelihood_calls++;
    }
}

/** Copy an array of points into preallocated vectors for the x and y values */
static inline void
blit_points_to_arrays(const point_t points[], const size_t n,
                      double *t, double *l)
{
    size_t i;
    for(i = 0; i < n; i++, points++) {
        *t++ = points->t;
        *l++ = points->ll;
    }
}

double
estimate_ml_t(log_like_function_t *log_like, double t[],
              const size_t n_pts, const double tolerance, bsm_t* model,
              bool* success)
{
    *success = false;

    point_t *starting_pts = malloc(sizeof(point_t) * n_pts);
    evaluate_ll(log_like, t, n_pts, starting_pts);

    size_t n = n_pts;
    point_t* points = select_points(log_like, starting_pts, &n,
                                    DEFAULT_MAX_POINTS);
    free(starting_pts);

    if (points == NULL) {
        fprintf(stderr, "ERROR: select_points returned NULL\n");
        *success = false;
        return NAN;
    }

    curve_type_t m = classify_curve(points, n_pts);

    if (m == CRV_MONO_DEC) {
        const double ml_t = points[0].t;
        free(points);
        *success = true;
        return ml_t;
    } else if (m != CRV_ENC_MAXIMA) {
        fprintf(stderr, "ERROR: selected points don't enclose a maximum\n");
        free(points);
        *success = false;
        return NAN;
    }

    assert(n >= n_pts);
    if (n > n_pts) {
        /* Subset to top n_pts */
        subset_points(points, n, n_pts);
    }


#ifdef VERBOSE
    fprintf(stderr, "starting iterative fit\n");
#endif /* VERBOSE */

    /* Allocate an extra point for scratch */
    points = realloc(points, sizeof(point_t) * (n_pts + 1));

    size_t iter = 0;
    const point_t* max_pt = NULL;
    double* l = malloc(sizeof(double) * n_pts);
    double ml_t = 0.0;

    for(iter = 0; iter < MAX_ITERS; iter++) {
        max_pt = max_point(points, n_pts);

        /* Re-fit */
        lcfit_bsm_rescale(max_pt->t, max_pt->ll, model);
        blit_points_to_arrays(points, n_pts, t, l);
        lcfit_fit_bsm(n_pts, t, l, model, 250);

        ml_t = lcfit_bsm_ml_t(model);

        if(isnan(ml_t)) {
            fprintf(stderr,
                    "ERROR: lcfit_bsm_ml_t returned NaN"
                    ", model = { %.3f, %.3f, %.6f, %.6f }\n",
                    model->c, model->m, model->r, model->b);
            *success = false;
            break;
        }

        /* If the BSM ml_t is zero or negative, add a new, smaller
         * sample point instead. */
        if(ml_t == 0.0) {
            ml_t = points[0].t / 10.0;
        }

        /* convergence check */
        if(fabs(ml_t - max_pt->t) <= tolerance) {
            *success = true;
            break;
        }

        /* Warn if ml_t is outside the bracketed window. */
        size_t max_idx = max_pt - points;
        if(ml_t < points[max_idx - 1].t || ml_t > points[max_idx + 1].t) {
            fprintf(stderr,
                    "WARNING: BSM ml_t (%g) is outside bracketed window [%g, %g]"
                    ", model = { %.3f, %.3f, %.6f, %.6f }\n",
                    ml_t, points[max_idx - 1].t, points[max_idx + 1].t,
                    model->c, model->m, model->r, model->b);
        }

        points[n_pts].t = ml_t;
        points[n_pts].ll = log_like->fn(ml_t, log_like->args);
        ml_likelihood_calls++;

        /* Retain top n_pts by log-likelihood */
        sort_by_t(points, n_pts + 1);

        if(classify_curve(points, n_pts + 1) != CRV_ENC_MAXIMA) {
            fprintf(stderr, "ERROR: after iteration points no longer enclose a maximum\n");
            *success = false;
            ml_t = NAN;
            break;
        }
        subset_points(points, n_pts + 1, n_pts);
    }

    if(iter == MAX_ITERS) {
        fprintf(stderr, "WARNING: maximum number of iterations reached\n");
        *success = false;
    }

    free(l);
    free(points);

#ifdef VERBOSE
    fprintf(stderr, "ending iterative fit after %zu iteration(s)\n", iter+1);
#endif /* VERBOSE */
    return ml_t;
}
