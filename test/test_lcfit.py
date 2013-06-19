#!/usr/bin/env python
from ctypes import c_double, c_ulong, POINTER, byref, Structure, cdll, \
        CFUNCTYPE, c_bool
import math
import unittest

liblcfit = cdll.LoadLibrary("_build/debug/lcfit_src/liblcfit.so")


# Definitions
class bsm_t(Structure):
    _fields_ = [("c", c_double),
                ("m", c_double),
                ("r", c_double),
                ("b", c_double)]


class point_t(Structure):
    _fields_ = [("t", c_double),
                ("ll", c_double)]


class log_like_func_t(Structure):
    LL_FUNC = CFUNCTYPE(c_double, c_double, POINTER(None))
    _fields_ = [("fn", LL_FUNC),
                ("args", POINTER(None))]


CRV_UNKNOWN, CRV_MONO_INC, CRV_MONO_DEC, CRV_ENC_MINIMA, CRV_ENC_MAXIMA = range(5)

class ClassifyCurveTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 4
        self.point_array = point_t * self.n

    def test_decreasing(self):
        p = self.point_array((0.0, 1.0), (1.0, 0.8), (2.0, 0.3), (3.0, 0.0))
        result = liblcfit.classify_curve(p, self.n)
        self.assertEqual(result, CRV_MONO_DEC)

    def test_decreasing2(self):
        p = self.point_array((0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (3.0, 0.0))
        result = liblcfit.classify_curve(p, self.n)
        self.assertEqual(result, CRV_MONO_DEC)

    def test_increasing(self):
        p = self.point_array((0.0, 1.0), (1.0, 10), (2.0, 100), (3.0, 1000.0))
        result = liblcfit.classify_curve(p, self.n)
        self.assertEqual(result, CRV_MONO_INC)

    def test_enc_maxima(self):
        p = self.point_array((0.0, 1.0), (1.0, 10), (2.0, 100), (3.0, 10))
        result = liblcfit.classify_curve(p, self.n)
        self.assertEqual(result, CRV_ENC_MAXIMA)

    def test_enc_minima1(self):
        p = self.point_array((0.0, -1.0), (1.0, -10), (2.0, -100), (3.0, 10))
        result = liblcfit.classify_curve(p, self.n)
        self.assertEqual(result, CRV_ENC_MINIMA)

    def test_enc_minima2(self):
        n = 8
        point_array = point_t * n
        p = point_array((1.0000000000000001e-05, -33047.506890500852),
                        (0.0001, -33047.506890500852),
                        (0.001, -33047.506890500852),
                        (0.01, -33047.506890500852),
                        (0.10000000000000001, -33047.506890500867),
                        (0.30000000000000004, -33047.50689050091),
                        (0.5, -33047.506890500918),
                        (1, -33047.506890500888))
        result = liblcfit.classify_curve(p, n)
        self.assertEqual(result, CRV_ENC_MINIMA)


DEFAULT_MODEL = bsm_t(1200, 300, 1, 0.2)


def py_ll(t, v=None):
    """
    ML branch length is 2.1978
    """
    md = DEFAULT_MODEL
    c = md.c
    m = md.m
    r = md.r
    b = md.b

    expterm = math.exp(-r * (t + b))
    return c * math.log((1 + expterm) / 2.) + m * math.log((1 - expterm) * 2.)


class ChoosePointsTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ll = log_like_func_t(log_like_func_t.LL_FUNC(py_ll), None)
        cls.select_points = liblcfit.select_points
        cls.select_points.restype = POINTER(point_t)
        cls.select_points.argtypes = [POINTER(log_like_func_t),
                                      POINTER(point_t),
                                      POINTER(c_ulong),
                                      c_ulong]

    def test_basic(self):
        n = 3
        args = [(t, py_ll(t)) for t in (0.5, 1.0, 1.1)]
        init_array = (point_t * n)(*args)

        n_pts = c_ulong(n)
        max_pts = c_ulong(8)
        points = self.select_points(byref(self.ll), init_array,
                                    byref(n_pts), max_pts)

        self.assertEqual(4, n_pts.value)
        expected_x = [0.05, 0.5, 1.0, 1.1]
        for i in xrange(n_pts.value):
            self.assertAlmostEqual(expected_x[i], points[i].t)
            self.assertAlmostEqual(py_ll(expected_x[i]), points[i].ll)
        liblcfit.free(points)


class SortPointsTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Define functions
        cls.sort_by_like = liblcfit.sort_by_like
        cls.sort_by_like.restype = None
        cls.sort_by_like.argtypes = [POINTER(point_t), c_ulong]

        cls.sort_by_t = liblcfit.sort_by_t
        cls.sort_by_t.restype = None
        cls.sort_by_t.argtypes = [POINTER(point_t), c_ulong]

    def setUp(self):
        self.n = 4
        self.array_t = point_t * self.n

    def test_sort_like(self):
        arr = self.array_t((1, -4), (0.5, -500), (0.75, -0.75), (2.0, -28))
        self.sort_by_like(arr, self.n)

        t = [(float(i.t), float(i.ll)) for i in arr]
        self.assertEqual(sorted(t, key=lambda x: x[1], reverse=True), t)

    def test_sort_t(self):
        arr = self.array_t((1, -4), (0.5, -500), (0.75, -0.75), (2.0, -28))
        self.sort_by_t(arr, self.n)

        t = [(float(i.t), float(i.ll)) for i in arr]
        self.assertEqual(sorted(t, key=lambda x: x[0]), t)


class EstimateMLTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ll = log_like_func_t(log_like_func_t.LL_FUNC(py_ll), None)
        cls.estimate_ml_t = liblcfit.estimate_ml_t
        cls.estimate_ml_t.restype = c_double
        cls.estimate_ml_t.argtypes = [POINTER(log_like_func_t),
                                      POINTER(c_double),
                                      c_ulong,
                                      c_double,
                                      POINTER(bsm_t),
                                      POINTER(c_bool)]
        cls.lcfit_bsm_ml_t = liblcfit.lcfit_bsm_ml_t
        cls.lcfit_bsm_ml_t.restype = c_double
        cls.lcfit_bsm_ml_t.argtypes = [POINTER(bsm_t)]

    def test_basic(self):
        n = 4
        bsm = bsm_t(1800.0, 400.0, 1.0, 0.5)
        pts = (c_double * n)(0.1, 0.5, 1.0, 1.5)
        success = c_bool()
        r = self.estimate_ml_t(self.ll, pts, n, 1e-3, byref(bsm),
                               byref(success))
        self.assertEqual(True, bool(success))
        self.assertAlmostEqual(self.lcfit_bsm_ml_t(byref(DEFAULT_MODEL)), r,
                               places=2)

    def test_no_max_enclosed(self):
        n = 4
        bsm = bsm_t(1800.0, 400.0, 1.0, 0.5)
        pts = (c_double * n)(1.0, 1.1, 1.4, 1.5)
        success = c_bool()
        r = self.estimate_ml_t(self.ll, pts, n, 1e-3, byref(bsm),
                               byref(success))
        self.assertEqual(True, bool(success))
        self.assertAlmostEqual(self.lcfit_bsm_ml_t(byref(bsm)), r,
                               places=2)


class SubsetPointsTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.subset_points = liblcfit.subset_points
        cls.subset_points.restype = None
        cls.subset_points.argtypes = [POINTER(point_t), c_ulong, c_ulong]

    def _test(self, expected, points, k):
        n = len(points)
        assert len(expected) >= k
        arr = (point_t * n)(*points)
        self.subset_points(arr, n, k)
        res = [(arr[i].t, arr[i].ll) for i in xrange(n)]
        self.assertEqual(expected, res)

    def test_n_equals_k(self):
        pts = [(0.1, 0.4), (0.2, 0.5), (0.3, 0.6), (0.4, 0.1)]
        self._test(pts, pts, len(pts))

    def test_n_lt_k1(self):
        p = [(0.1, 0.4), (0.2, 0.5), (0.3, 0.6), (0.4, 0.1)]
        expected_order = [1, 2, 3, 0]
        expected = [p[i] for i in expected_order]
        self._test(expected, p, 3)

    def test_n_lt_k2(self):
        p = [(0.1, 0.4), (0.2, 0.5), (0.3, 0.3), (0.4, 0.1)]
        expected_order = [0, 1, 2, 3]
        expected = [p[i] for i in expected_order]
        self._test(expected, p, 3)

    def test_n_lt_k3(self):
        p = [(0.1, 0.4), (0.2, 0.5), (0.3, 0.6), (0.4, 0.8),
             (0.45, 0.74), (0.5, 0.6)]
        expected_order = [2, 3, 4, 5, 1, 0]
        expected = [p[i] for i in expected_order]
        self._test(expected, p, 4)

if __name__ == '__main__':
    unittest.main()
