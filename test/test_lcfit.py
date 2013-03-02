#!/usr/bin/env python
import ctypes
import math
import unittest

liblcfit = ctypes.cdll.LoadLibrary("_build/debug/liblcfit-shared.so")


class point_t(ctypes.Structure):
    _fields_ = [("x", ctypes.c_double),
                ("y", ctypes.c_double)]

MONO_UNKNOWN = 0
MONO_INC = 1
MONO_DEC = 2
NON_MONOTONIC = 3


class MonotonicityTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 4
        self.point_array = point_t * self.n

    def test_decreasing(self):
        p = self.point_array((0.0, 1.0), (1.0, 0.8), (2.0, 0.3), (3.0, 0.0))
        result = liblcfit.monotonicity(p, self.n)
        self.assertEqual(result, MONO_DEC)

    def test_increasing(self):
        p = self.point_array((0.0, 1.0), (1.0, 10), (2.0, 100), (3.0, 1000.0))
        result = liblcfit.monotonicity(p, self.n)
        self.assertEqual(result, MONO_INC)

    def test_non_monotonic(self):
        p = self.point_array((0.0, 1.0), (1.0, 10), (2.0, 100), (3.0, 10))
        result = liblcfit.monotonicity(p, self.n)
        self.assertEqual(result, NON_MONOTONIC)

# Machinery for log-likelihood
LL_FUNC = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.POINTER(None))


def py_ll(t, v=None):
    c = 1200.0
    m = 1000.0
    r = 1.0
    b = 0.2

    expterm = math.exp(-r * (t + b))
    return c * math.log((1 + expterm) / 2.) + m * math.log((1 - expterm) * 2.)


ll = LL_FUNC(py_ll)


class ChoosePointsTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.select_points = liblcfit.select_points
        cls.select_points.restype = ctypes.POINTER(point_t)
        cls.select_points.argtypes = [LL_FUNC, ctypes.POINTER(ctypes.c_double),
                                      ctypes.POINTER(ctypes.c_ulong), ctypes.c_ulong,
                                      ctypes.POINTER(None)]

    def test_basic(self):
        n = 3
        init_array_t = ctypes.c_double * n
        init_array = init_array_t(0.5, 1.0, 1.1)

        n = ctypes.c_ulong(n)
        max_pts = ctypes.c_ulong(8)
        points = self.select_points(ll, init_array, ctypes.byref(n), max_pts, None)

        self.assertEqual(5, n.value)
        expected_x = [0.5, 1.0, 1.1, 2.2, 4.4]
        for i in xrange(n.value):
            self.assertAlmostEqual(expected_x[i], points[i].x)
            self.assertAlmostEqual(py_ll(expected_x[i]), points[i].y)


class SortPointsTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sort_by_like = liblcfit.sort_by_like
        cls.sort_by_like.restype = None
        cls.sort_by_like.argtypes = [ctypes.POINTER(point_t), ctypes.c_ulong]

    def test_sort(self):
        n = 4
        array_t = point_t * n
        arr = array_t((1, -4), (0.5, -500), (0.75, -0.75), (2.0, -28))
        self.sort_by_like(arr, n)

        t = [(float(i.x), float(i.y)) for i in arr]
        self.assertEqual(sorted(t, key=lambda x: x[1], reverse=True), t)

if __name__ == '__main__':
    unittest.main()
