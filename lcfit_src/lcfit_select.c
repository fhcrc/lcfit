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

#include "lcfit2.h"

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

static void print_points(FILE* fp, const point_t* points, const size_t n)
{
    char* sep = "";
    for (size_t i = 0; i < n; ++i) {
        fprintf(fp, "%s%g", sep, points[i].t);
        sep = ", ";
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

double bound_point(const double proposed_t, const point_t* points,
                   const size_t n_pts, const double min_t, const double max_t)
{
    double next_t = proposed_t;

    // Ensure the next branch length to evaluate is within bounds.
    if (next_t < min_t) {
        next_t = min_t;
    } else if (next_t > max_t) {
        next_t = max_t;
    }

    // If the next branch length is equal to the minimum or
    // maximum of the already-evaluated branch lengths, split the
    // difference between it and its neighbor instead.
    if (next_t == points[0].t) {
        next_t = points[0].t + (points[1].t - points[0].t) / 2.0;
    } else if (next_t == points[n_pts - 1].t) {
        next_t = points[n_pts - 2].t +
                 (points[n_pts - 1].t - points[n_pts - 2].t) / 2.0;
    }

    return next_t;
}

point_t*
select_points(log_like_function_t *log_like, const point_t starting_pts[],
              size_t *num_pts, const size_t max_pts, const double min_t,
              const double max_t)
{
    size_t n = *num_pts;
    assert(n >= 3);

    point_t* points = malloc(sizeof(point_t) * max_pts);

    /* Copy already-evaluated points */
    memcpy(points, starting_pts, sizeof(point_t) * n);

    /* Add additional samples until the evaluated branch lengths enclose a
     * maximum or the maximum number of points is reached. */
    for (; n < max_pts; ++n) {
        curve_type_t curvature = classify_curve(points, n);

        if (curvature == CRV_ENC_MAXIMA) {
            break;
        } else if (curvature == CRV_ENC_MINIMA || curvature == CRV_UNKNOWN) {
            free(points);
            return NULL;
        }

        double proposed_t = 0.0;

        if (curvature == CRV_MONO_INC) {
            proposed_t = points[n - 1].t * 2.0;
        } else { /* curvature == CRV_MONO_DEC */
            proposed_t = points[0].t / 10.0;
        }

        double next_t = bound_point(proposed_t, points, n, min_t, max_t);

        points[n].t = next_t;
        points[n].ll = log_like->fn(next_t, log_like->args);
        bracket_likelihood_calls++;

        sort_by_t(points, n + 1);
    }

    *num_pts = n;

    if (n < max_pts) {
        return realloc(points, sizeof(point_t) * n);
    }

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
    sort_by_t(p, n);

    if (k == n) return;
    assert(k < n);

    curve_type_t curvature = classify_curve(p, n);
    assert(curvature == CRV_ENC_MAXIMA || curvature == CRV_MONO_DEC);

    if (curvature == CRV_MONO_DEC) {
        /* Keep leftmost k points, which are already in place. */
        return;
    }

    /* Otherwise, the curvature is CRV_ENC_MAXIMA. Keep the maximum
     * and its neighbors and sort the rest by likelihood. */
    assert(n >= 3);

    size_t max_idx = point_max_index(p, n);

    if (max_idx > 1) {
        const size_t n_before = max_idx - 1;
        const size_t s = n_before * sizeof(point_t);
        point_t *buf = malloc(s);
        memcpy(buf, p, s);
        memmove(p, p + max_idx - 1, sizeof(point_t) * 3);
        memcpy(p + 3, buf, s);
        free(buf);
    }

    if (k > 3) {
        sort_by_like(p + 3, n - 3);
    }

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

double rel_err(double expected, double actual)
{
    return fabs((expected - actual) / expected);
}

double
estimate_ml_t(log_like_function_t *log_like, const double* t,
              size_t n_pts, const double tolerance, bsm_t* model,
              bool* success, const double min_t, const double max_t)
{
    *success = false;

    point_t *starting_pts = malloc(sizeof(point_t) * n_pts);
    evaluate_ll(log_like, t, n_pts, starting_pts);

    const size_t orig_n_pts = n_pts;
    point_t* points = select_points(log_like, starting_pts, &n_pts,
                                    DEFAULT_MAX_POINTS, min_t, max_t);
    free(starting_pts);

    if (points == NULL) {
        fprintf(stderr, "ERROR: select_points returned NULL\n");
        *success = false;
        return NAN;
    }

    curve_type_t curvature = classify_curve(points, n_pts);

    if (!(curvature == CRV_ENC_MAXIMA || curvature == CRV_MONO_DEC)) {
        fprintf(stderr, "ERROR: "
                "points don't enclose a maximum and aren't decreasing\n");

        free(points);
        *success = false;
        return NAN;
    }

    /* From here on, curvature is CRV_ENC_MAXIMA or CRV_MONO_DEC, and
     * thus ml_t is zero or positive (but not infinite). */

    assert(n_pts >= orig_n_pts);
    if (n_pts > orig_n_pts) {
        /* Subset to top orig_n_pts */
        subset_points(points, n_pts, orig_n_pts);
        n_pts = orig_n_pts;
    }

    assert(points[0].t >= min_t);
    assert(points[n_pts - 1].t <= max_t);

#ifdef VERBOSE
    fprintf(stderr, "starting iterative fit\n");

    fprintf(stderr, "starting points: ");
    print_points(stderr, points, n_pts);
    fprintf(stderr, "\n");
#endif /* VERBOSE */

    /* Allocate an extra point for scratch */
    points = realloc(points, sizeof(point_t) * (n_pts + 1));

    size_t iter = 0;
    const point_t* max_pt = NULL;
    double ml_t = 0.0;
    double prev_t = 0.0;

    for (iter = 0; iter < MAX_ITERS; iter++) {
        max_pt = max_point(points, n_pts);

        /* Re-fit */
        lcfit_bsm_rescale(max_pt->t, max_pt->ll, model);

        double* tbuf = malloc(sizeof(double) * n_pts);
        double* lbuf = malloc(sizeof(double) * n_pts);

        blit_points_to_arrays(points, n_pts, tbuf, lbuf);

        double* wbuf = malloc(sizeof(double) * n_pts);
        double alpha = (double) iter / (MAX_ITERS - 1);

        for (size_t i = 0; i < n_pts; ++i) {
            wbuf[i] = exp(alpha * (lbuf[i] - max_pt->ll));
        }

#ifdef VERBOSE
        fprintf(stderr, "weights: ");
        for (size_t i = 0; i < n_pts; ++i) {
            fprintf(stderr, "%g ", wbuf[i]);
        }
        fprintf(stderr, "\n");
#endif /* VERBOSE */

        lcfit_fit_bsm_weight(n_pts, tbuf, lbuf, wbuf, model, 250);

        free(tbuf);
        free(lbuf);
        free(wbuf);

        ml_t = lcfit_bsm_ml_t(model);

        if (isnan(ml_t)) {
            fprintf(stderr, "ERROR: "
                    "lcfit_bsm_ml_t returned NaN"
                    ", model = { %.3f, %.3f, %.6f, %.6f }\n",
                    model->c, model->m, model->r, model->b);
            *success = false;
            break;
        }

        if (curvature == CRV_ENC_MAXIMA) {
            /* Stop if the modeled maximum likelihood branch length is
             * within tolerance of the empirical maximum. */
            if (rel_err(max_pt->t, ml_t) <= tolerance) {
                *success = true;
                break;
            }
        }

        double next_t = bound_point(ml_t, points, n_pts, min_t, max_t);

        /* Stop if the next sample point is within tolerance of the
         * previous sample point. */
        if (rel_err(prev_t, next_t) <= tolerance) {
            *success = true;
            break;
        }

        points[n_pts].t = next_t;
        points[n_pts].ll = log_like->fn(next_t, log_like->args);
        ml_likelihood_calls++;

        prev_t = next_t;

        sort_by_t(points, n_pts + 1);
        curvature = classify_curve(points, n_pts + 1);

        if (!(curvature == CRV_ENC_MAXIMA || curvature == CRV_MONO_DEC)) {
            fprintf(stderr, "ERROR: "
                    "after iteration points don't enclose a maximum "
                    "and aren't decreasing\n");
            *success = false;
            break;
        }

        ++n_pts;
        points = realloc(points, sizeof(point_t) * (n_pts + 1));

#ifdef VERBOSE
        fprintf(stderr, "current points: ");
        print_points(stderr, points, n_pts);
        fprintf(stderr, "\n");
#endif /* VERBOSE */
    }

    if (iter == MAX_ITERS) {
        fprintf(stderr, "WARNING: maximum number of iterations reached\n");
    }

    free(points);

#ifdef VERBOSE
    fprintf(stderr, "ending iterative fit after %zu iteration(s)\n", iter);
#endif /* VERBOSE */

    if (ml_t < min_t) {
        ml_t = min_t;
    } else if (ml_t > max_t) {
        ml_t = max_t;
    }

    return ml_t;
}

double lcfit_fit_auto(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                      bsm_t* model, const double min_t, const double max_t)
{
    double d1;
    double d2;
    double t0 = lcfit_maximize(lnl_fn, lnl_fn_args, min_t, max_t, &d1, &d2);

    if (fabs(d1) < 0.1) {
        lcfit2_bsm_t lcfit2_model = {1100.0, 800.0, t0, d1, d2};
        const double alpha = 0.0;

        lcfit2_fit_auto(lnl_fn, lnl_fn_args, &lcfit2_model, min_t, max_t, alpha);
        lcfit2_to_lcfit4(&lcfit2_model, model);
    } else {
        // HACK: basically copied from AdHocIntegrator.cpp

        log_like_function_t lnl_fn_wrapper = {lnl_fn, lnl_fn_args};
        double t[4] = {0.1, 0.15, 0.5, 1.0};
        const double tolerance = 1e-3;
        bool success = false;

        t0 = estimate_ml_t(&lnl_fn_wrapper, t, 4, tolerance, model, &success, min_t, max_t);
    }

    return t0;
}
