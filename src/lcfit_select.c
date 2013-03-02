#include "lcfit_select.h"
#include <assert.h>
#include <string.h>

#define FALSE 0
#define TRUE 1

monotonicity_t
monotonicity(const point_t* points, const size_t n)
{
    short maybe_inc = 1, maybe_dec = 1;
    size_t i;
    const point_t *last = points;
    for(i = 1; i < n; ++i) {
        const point_t *p = &points[i];

        assert(p->x >= last->x);
        if(p->y > last->y)
            maybe_dec = 0;
        else if(p->y < last->y)
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
select_points(log_like_function_t log_like, const double* starting_pts, size_t *num_pts, const size_t max_pts, void* args)
{
    point_t* points = malloc(sizeof(point_t) * (max_pts));
    size_t i, n = *num_pts;

    for(i = 0; i < max_pts; ++i) {
        points[i].x = -1; points[i].y = -1;
    }
    for(i = 0; i < n; ++i) {
        points[i].x = starting_pts[i];
        points[i].y = log_like(points[i].x, args);
    }

    /* Add additional samples until the evaluated branch lengths enclose a
     * maximum. */
    size_t offset = 0;  /* Position */
    double d = 0.0;     /* Branch length */

    monotonicity_t c = monotonicity(points, n);
    do {
        switch(c) {
            case NON_MONOTONIC:
                /* Add a point between the first and second try */
                d = (points[1].x + points[2].x) / 2.0;
                offset = 1;
                /* shift */
                memmove(points + offset, points + offset + 1, sizeof(point_t) * (n - offset));
                break;
            case MONO_INC:
                /* Double largest value */
                d = points[n - 1].x * 2.0;
                offset = n;
                break;
            case MONO_DEC:
                /* Add new smallest value - order of magnitude lower */
                d = points[0].x / 10.0;
                offset = 0;
                /* shift */
                memmove(points, points + 1, sizeof(point_t) * n);
                break;
            default:
                assert(FALSE);
        }

        const double l = log_like(d, args);
        points[offset].x = d;
        points[offset].y = l;
        n++;

        c = monotonicity(points, n);
    } while(n <= max_pts && c != NON_MONOTONIC);

    *num_pts = n;

    if(n < max_pts)
        return realloc(points, sizeof(point_t) * n);
    else
        return points;
}

int dec_like_cmp(const void *p1, const void *p2)
{
    const point_t *pt1 = (const point_t*) p1,
                  *pt2 = (const point_t*) p2;
    if(pt1->y < pt2->y) return 1;
    if(pt1->y == pt2->y) return 0;
    return -1;
}

void sort_by_like(point_t points[], const size_t n)
{
    qsort(points, n, sizeof(point_t), dec_like_cmp);
}
