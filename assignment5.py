"""
In this assignment you should fit a model function of your choice to data 
that you sample from a contour of given shape. Then you should calculate
the area of that shape. 

The sampled data is very noisy so you should minimize the mean least squares 
between the model you fit and the data points you sample.  

During the testing of this assignment running time will be constrained. You
receive the maximal running time as an argument for the fitting method. You 
must make sure that the fitting function returns at most 5 seconds after the 
allowed running time elapses. If you know that your iterations may take more 
than 1-2 seconds break out of any optimization loops you have ahead of time.

Note: You are allowed to use any numeric optimization libraries and tools you want
for solving this assignment. 
Note: !!!Despite previous note, using reflection to check for the parameters 
of the sampled function is considered cheating!!! You are only allowed to 
get (x,y) points from the given shape by calling sample(). 
"""
import math
import operator
import warnings
from scipy import interpolate
from scipy.interpolate import BSpline, splev, splprep
from functionUtils import AbstractShape
from functools import cmp_to_key, reduce
import numpy as np


class MyShape(AbstractShape):

    # change this class with anything you need to implement the shape

    def __init__(self, area: float, func: callable):

        self._area = area
        self._func = func
        self.ass = Assignment5()

    def area(self):

        if self._area is None:
            return self.ass.area(self._func)
        return self._area

    def contour(self, n):
        return self._func(n)
        pass


class Assignment5:
    def __init__(self):
        """
        Here goes any one time calculation that need to be made before
        solving the assignment for specific functions.
        """

        pass
    def interp1(self, points,n):
        coords = points.copy()

        center = tuple(
            map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), coords), [len(coords)] * 2))
        pts = (sorted(coords, key=lambda coord: (-135 - math.degrees(
            math.atan2(*tuple(map(operator.sub, coord, center))[::-1]))) % 360))

        pad = 3
        pts = np.pad(pts, [(pad, pad), (0, 0)], mode='wrap')
        x, y = pts.T
        i = np.arange(0, len(pts))

        interp_i = np.linspace(pad, i.max() - pad + 1, 5 * (i.size - 2 * pad))
        xi = interpolate.interp1d(i, x, kind='cubic')(interp_i)
        yi = interpolate.interp1d(i, y, kind='cubic')(interp_i)
        return [(x, y) for x, y in zip(xi, yi)]

    def interp2(self, points, n):
        x_points, y_points = [], []
        for pt in points:
            x_points.append(pt[0])
            y_points.append(pt[1])
        x_points.sort()
        y_points.sort()
        t, c, k = interpolate.splrep(x_points, y_points, s=0, k=3)
        spline = BSpline(t, c, k)
        X = np.linspace(0, 1, num=n)
        ans = [np.array((x, spline(x)),  dtype=np.longdouble) for x in X]
        return ans


    def interp5(self,points,n):
        x = [(p[0],p[1]) for p in points]
        x = np.asarray(x)
        xs = (x - x.mean(0))
        x_sort = xs[np.angle((xs[:, 0] + 1j * xs[:, 1])).argsort()]
        tck, u = splprep(x_sort.T, u=None, s=0.0, per=1)
        u_new = np.linspace(u.min(), u.max(), 1000)
        xi, yi = splev(u_new, tck, der=0)
        return [(x,y) for x,y in zip(xi,yi)]


    def get_center(self, points):
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        cen_x = sum(x) / len(points)
        cen_y = sum(y) / len(points)
        left_down, top_left, top_right, down_right = [], [], [], []
        for x, y in points:
            if x > cen_x and y < cen_y:
                down_right.append((x, y))
            elif x > cen_x and y > cen_y:
                top_right.append((x, y))
            elif x < cen_x and y < cen_y:
                left_down.append((x, y))
            else:
                top_left.append((x, y))
            left_down = sorted(left_down, key=cmp_to_key(lambda p1, p2: p1[0] - p2[0]))
            down_right = sorted(down_right, key=cmp_to_key(lambda p1, p2: p1[0] - p2[0]))
            top_right = sorted(top_right, key=cmp_to_key(lambda p1, p2: p1[0] - p2[0]))
            top_left = sorted(top_left, key=cmp_to_key(lambda p1, p2: p1[0] - p2[0]))
        return left_down, top_left, top_right, down_right

    def area(self, contour: callable, maxerr=0.001) -> np.float32:
        """
        Compute the result of the shape with the given contour.

        Parameters
        ----------
        contour : callable
            Same as AbstractShape.contour
        maxerr : TYPE, optional
            The target error of the result computation. The default is 0.001.

        Returns
        -------
        The result of the shape.

        """

        n = 10000
        pts = contour(n)
        n = len(pts)
        result = np.float128(0)
        if n<2:
            return np.float32(0)
        prev = n - 1
        for i in range(0, n):
            result += np.float128((pts[prev][0] + pts[i][0]) * (pts[prev][1] - pts[i][1]))
            prev = i
        result = np.float128(result / 2)
        return np.abs(np.float32(result))


    def fit_shape(self, sample: callable, maxtime: float) -> AbstractShape:

        """
        Build a function that accurately fits the noisy data points sampled from
        some closed shape.

        Parameters
        ----------
        sample : callable.
            An iterable which returns a data point that is near the shape contour.
        maxtime : float
            This function returns after at most maxtime seconds.

        Returns
        -------
        An object extending AbstractShape.
        """
        warnings.filterwarnings('ignore')

        Tt = time.time()
        def creation(left_down, down_right, top_right, top_left, pts):
            points = []
            points.extend(top_left)
            points.extend(top_right)
            points.extend(down_right)
            points.extend(left_down)
            def func_contour(n):
                if len(points) < 2:
                    return []
                len1 = len(left_down)
                len2 = len(down_right)
                len3 = len(top_right)
                len4 = len(top_left)
                mid = ((len1 + len2 + len3 + len4) / 4)
                if mid > len1 * 10 or mid > len2 * 10 or mid > len3 * 10 or mid > len4 * 10:
                    return self.interp2(points, n)
                return self.interp5(points,n)
            return func_contour



        # replace these lines with your solution
        T = time.time()
        coords = []
        sum_x = 0
        sum_y = 0
        n = 0
        while n < 30000:
            if time.time() - Tt > maxtime - 2:
                break
            value = sample()
            n += 1
            coords.append(value)
            sum_x = sum_x + value[0]
            sum_y = sum_y + value[1]
        if n != 0:
            cen_x = sum_x / n
            cen_y = sum_y / n
        else:
            cen_x = 0
            cen_y = 0

        left_down, top_left, top_right, down_right = [], [], [], []
        for x, y in coords:
            if x > cen_x and y < cen_y:
                down_right.append((x, y))
            elif x > cen_x and y > cen_y:
                top_right.append((x, y))
            elif x < cen_x and y < cen_y:
                left_down.append((x, y))
            else:
                top_left.append((x, y))

        """"
        print(len(left_down))
        print(len(top_left))
        print(len(top_right))
        print(len(down_right))
        if len(top_right) < 2000:
            left_down, top_left, top_right, down_right = self.get_center(left_down)
            left_down, top_left, top_right, down_right = self.get_center(top_left)
            left_down, top_left, top_right, down_right = self.get_center(down_right)
            _func1 = creation(left_down, down_right, top_right, top_left)
            _func2 = creation(left_down, down_right, top_right, top_left)
            _func3 = creation(left_down, down_right, top_right, top_left)
            shape_area = self.area(_func1) + self.area(_func2) + self.area(_func3)
            return MyShape(shape_area, _func1)
            """

        shape_area = None
        if time.time() - T >= maxtime - 3:
            func_c = creation(left_down, down_right, top_right, top_left,coords)
            return MyShape(shape_area, func_c)

        _func = creation(left_down, down_right, top_right, top_left,coords)

        shape_area = abs(self.area(_func))
        return MyShape(shape_area, _func)
        # replace these lines with your solution


##########################################################################


from sampleFunctions import *


class TestAssignment4(unittest.TestCase):
    """
    def test_return(self):
        circ = noisy_circle(cx=1, cy=1, radius=1, noise=0.1)
        ass4 = Assignment5()
        T = time.time()
        shape = ass4.fit_shape(sample=circ, maxtime=5)
        T = time.time() - T
        self.assertTrue(isinstance(shape, AbstractShape))
        self.assertLessEqual(T, 5)
        """

    def test_circle_area(self):
        circ = noisy_circle(cx=1, cy=1, radius=1, noise=0.1)
        ass4 = Assignment5()
        T = time.time()
        shape = ass4.fit_shape(sample=circ, maxtime=5)
        T = time.time() - T
        a = shape.area()
        print("circle - ", abs(abs(a - np.pi)))
        self.assertLess(abs(a - np.pi), 0.01)
        self.assertLessEqual(T, 32)

    def test_bezier_fit(self):
        circ = noisy_circle(cx=1, cy=1, radius=1, noise=0.1)
        ass5 = Assignment5()
        T = time.time()
        shape = ass5.fit_shape(sample=circ, maxtime=3)
        T = time.time() - T
        a = shape.area()
        print(abs(a - np.pi))
        self.assertLess(abs(a - np.pi), 0.01)
        self.assertLessEqual(T, 32)

    def test_circle_area_from_contour(self):
        circ = Circle(cx=1, cy=1, radius=1, noise=0.0)
        ass4 = Assignment5()
        T = time.time()
        a_computed = ass4.area(contour=circ.contour, maxerr=0.1)
        T = time.time() - T
        a_true = circ.area()
        print(abs((a_true - a_computed)))
        self.assertLess(abs((a_true - a_computed) / a_true), 0.1)


if __name__ == "__main__":
    unittest.main()
