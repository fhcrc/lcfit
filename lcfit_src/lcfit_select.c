/**
 * \file lcfit_select.c
 * \brief Point selection implementation
 */
#include "lcfit_select.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define FALSE 0
#define TRUE 1

static const double DEFAULT_START[] = {0.1, 0.5, 1.0};
static const size_t N_DEFAULT_START = 3;
/** Default maximum number of points to evaluate in select_points */
static const size_t DEFAULT_MAX_POINTS = 8;

monotonicity_t
monotonicity(const point_t points[], const size_t n)
{
    const point_t *p = points;
    short maybe_inc = 1, maybe_dec = 1;
    size_t i;

    const point_t *last = p++;
    for(i = 1; i < n; ++i, ++p) {
        assert(p->t >= last->t && "Points not sorted!");
        if(p->ll > last->ll)
            maybe_dec = 0;
        else if(p->ll < last->ll)
            maybe_inc = 0;
        last = p;
    }

    if(!maybe_inc && !maybe_dec) return NON_MONOTONIC;
    else if(maybe_inc) return MONO_INC;
    else if(maybe_dec) return MONO_DEC;
    assert(FALSE && "Unknown monotonicity");
    return MONO_UNKNOWN;
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

    monotonicity_t c = monotonicity(points, n);
    do {
        switch(c) {
            case NON_MONOTONIC:
                /* Add a point between the first and second try */
                d = points[0].t + ((points[1].t - points[0].t) / 2.0);
                offset = 1;
                /* shift */
                memmove(points + offset + 1, points + offset,
                        sizeof(point_t) * (n - offset));
                break;
            case MONO_INC:
                /* Double largest value */
                d = points[n - 1].t * 2.0;
                offset = n;
                break;
            case MONO_DEC:
                /* Add new smallest value - order of magnitude lower */
                d = points[0].t / 10.0;
                offset = 0;
                /* shift */
                memmove(points + 1, points, sizeof(point_t) * n);
                break;
            default:
                assert(FALSE);
        }

        const double l = log_like->fn(d, log_like->args);
        points[offset].t = d;
        points[offset].ll = l;

        c = monotonicity(points, n++);
    } while(n < max_pts && c != NON_MONOTONIC);

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
    }
    return m;
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
              const size_t n_pts, const double tolerance, bsm_t* model)
{
    size_t iter = 0;
    const point_t *max_pt;
    double *l = malloc(sizeof(double) * n_pts);

    /* We allocate an extra point for scratch */
    point_t *points = malloc(sizeof(point_t) * (n_pts + 1));
    evaluate_ll(log_like, t, n_pts, points);

    /* First, classify points */
    monotonicity_t m = monotonicity(points, n_pts);

    /* If the function is no longer monotonic, start over */
    if(m != NON_MONOTONIC) {
        free(points);

        size_t n = N_DEFAULT_START;
        point_t *start_pts = malloc(sizeof(point_t) * n);
        evaluate_ll(log_like, DEFAULT_START, n, start_pts);
        points = select_points(log_like, start_pts, &n, DEFAULT_MAX_POINTS);
        free(start_pts);
        assert(n >= n_pts);
        if(n > n_pts) {
            /* Subset to top n_pts */
            sort_by_like(points, n);
            sort_by_t(points, n_pts);
        }


        /* Allocate an extra point for scratch */
        points = realloc(points, sizeof(point_t) * (n_pts + 1));
    }

    max_pt = max_point(points, n_pts);
    double ml_t = -DBL_MAX;

    /* Re-fit */
    lcfit_bsm_rescale(max_pt->t, max_pt->ll, model);
    blit_points_to_arrays(points, n_pts, t, l);
    lcfit_fit_bsm(n_pts, t, l, model);

    /* TODO: factor out magic number for max iters */
    for(iter = 0; iter < 100; iter++) {
        ml_t = lcfit_bsm_ml_t(model);

        /* convergence check */
        if(fabs(ml_t - max_pt->t) <= tolerance) {
            break;
        }

        /* Add ml_t estimate */
        if(ml_t < 0) {
          ml_t = 1e-8;
          break;
        }
        points[n_pts].t = ml_t;
        points[n_pts].ll = log_like->fn(ml_t, log_like->args);

        /* Retain top n_pts by log-likelihood */
        sort_by_like(points, n_pts + 1);
        sort_by_t(points, n_pts);

        blit_points_to_arrays(points, n_pts, t, l);
        lcfit_fit_bsm(n_pts, t, l, model);
    }

    free(l);
    free(points);

    return ml_t;
}
