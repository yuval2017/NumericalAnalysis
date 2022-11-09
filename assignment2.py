"""
In this assignment you should find the intersection points for two functions.
"""
import sys
from collections.abc import Iterable


def falsePosition(x0, x1, f, e):
    x2 = x0
    while abs(f(x2)) >= e:
        x2 = x0 - (x1 - x0) * f(x0) / (f(x1) - f(x0))
        if f(x0) * f(x2) < 0:
            x1 = x2
        else:
            x0 = x2
    return x2
def get_y(portion, f, maxerr=0.001):
    fx = None
    if f(portion[0]) * f(portion[1]) < 0:
        fx = falsePosition(portion[0], portion[1], f, maxerr)
    return fx


class Assignment2:
    def __init__(self):

        pass

    def intersections(self, f1: callable, f2: callable, a: float, b: float, maxerr=0.001) -> Iterable:
        """
        Find as many intersection points as you can. The assignment will be
        tested on functions that have at least two intersection points, one
        with a positive x and one with a negative x.
        
        This function may not work correctly if there is infinite number of
        intersection points. 


        Parameters
        ----------
        f1 : callable
            the first given function
        f2 : callable
            the second given function
        a : float
            beginning of the interpolation range.
        b : float
            end of the interpolation range.
        maxerr : float
            An upper bound on the difference between the
            function values at the approximate intersection points.


        Returns
        -------
        X : iterable of approximate intersection Xs such that for each x in X:
            |f1(x)-f2(x)|<=maxerr.

        """
        f = lambda a: f1(a) - f2(a)
        """"
        def intersections(a: float, b: float, maxerr=0.001) -> Iterable:
            eps = sys.float_info.epsilon
            ans = []
            p2 = a + maxerr
            div = 1
            found_root = False
            while a < b:
                if abs(f(a)) < maxerr:
                    found_root = True
                    ans.append(a)
                    curr_x = a + eps
                    div = 1
                else:
                    curr_x = get_y((a, p2), f, maxerr)
                    if (curr_x is not None):
                        found_root = True
                        div = 1
                        ans.append(curr_x)
                        curr_x = curr_x + eps
                    else:
                        found_root = False
                        curr_x = a
                        if div <= 2:
                            div = div + 1
                a = curr_x + maxerr / div
                p2 = a + maxerr
            if not found_root and abs(f(b)) < maxerr:
                ans.append(b)



            return ans
            """

        def creat_condition(a):
            if f(a) < 0:
                return lambda b: f(b) > 0
            else:
                return lambda b: f(b) < 0

        def next_root(left, b, maxrrr):
            right = left + maxrrr
            if (right > b):
                return None
            condition = creat_condition(left)
            while not condition(right):
                if abs(f(left)) < maxerr:
                    return left
                left = right
                right = right + (maxrrr)
                if (right > b):
                    return None
            return falsePosition(left, right, f, maxerr)

        def inter(a, b, maxerr):
            ans = []
            eps = sys.float_info.epsilon
            x = next_root(a, b, maxerr)
            while (x is not None):
                if x <= b:
                    ans.append(x)
                a = x + maxerr + eps
                x = next_root(a, b, maxerr)
            if len(ans) != 0 and abs(ans[len(ans) - 1] - b) < maxerr and abs(f(b)) < maxerr:
                ans.append(b)
                """
            for i in range(len(ans)):
                for j in range(len(ans)):
                    if (i != j) and (abs(ans[i] - ans[j]) <= maxerr):
                        print("False")
                        """
            return ans

        return inter(a, b, maxerr)


##########################################################################


from sampleFunctions import *


class TestAssignment2(unittest.TestCase):

    def test_sqr(self):

        ass2 = Assignment2()

        f1 = lambda x: np.sin(100 / x)
        f2 = lambda x: 0
        T = time.time()
        X = ass2.intersections(f1, f2, 1, 4, maxerr=0.001)
        T = time.time() - T
        print(T)
        print(X)
        print(len(X))
        for x in X:
            self.assertGreaterEqual(0.001, abs(f1(x) - f2(x)))

    def test_sqr2(self):
        print("test2")
        ass2 = Assignment2()
        f1, f2 = lambda x: x ** 2, lambda x: 0
        T = time.time()
        X = ass2.intersections(f1, f2, -1, 4, maxerr=0.001)
        T = time.time() - T
        print(T)
        print(X)
        print(len(X))
        for x in X:
            self.assertGreaterEqual(0.001, abs(f1(x) - f2(x)))

    def test_sqr1(self):
        print("test1")
        ass2 = Assignment2()
        f1, f2 = lambda x: -6 * (x - 7) * 4 + 24 * x - 29, lambda x: x * 2 - 10 * x - 5
        T = time.time()
        X = ass2.intersections(f1, f2, -100, 100, maxerr=0.001)
        T = time.time() - T
        print(T)
        print(X)
        print(len(X))
        for x in X:
            self.assertGreaterEqual(0.001, abs(f1(x) - f2(x)))

    def test_poly(self):

        ass2 = Assignment2()

        f1, f2 = randomIntersectingPolynomials(10)

        X = ass2.intersections(f1, f2, 1, 21, maxerr=0.001)
        print(X)

        for x in X:
            self.assertGreaterEqual(0.001, abs(f1(x) - f2(x)))


if __name__ == "__main__":
    unittest.main()
