"""
In this assignment you should find the area enclosed between the two given functions.
The rightmost and the leftmost x values for the integration are the rightmost and 
the leftmost intersection points of the two functions. 

The functions for the numeric answers are specified in MOODLE. 


This assignment is more complicated than Assignment1 and Assignment2 because: 
    1. You should work with float32 precision only (in all calculations) and minimize the floating point errors. 
    2. You have the freedom to choose how to calculate the area between the two functions. 
    3. The functions may intersect multiple times. Here is an example: 
        https://www.wolframalpha.com/input/?i=area+between+the+curves+y%3D1-2x%5E2%2Bx%5E3+and+y%3Dx
    4. Some of the functions are hard to integrate accurately. 
       You should explain why in one of the theoretical questions in MOODLE. 

"""
import math
import random
import sys
import numpy as np
from scipy.special import roots_legendre
import sampleFunctions


class Assignment3:
    def __init__(self):
        """
        Here goes any one time calculation that need to be made before 
        solving the assignment for specific functions. 
        """
    def gauss(self, a, b, f, n):
        [X,W] = roots_legendre(n+1)
        ans = np.float128(0)
        for x,w in zip(X,W):
            ans += np.float128(w * f(0.5 * (b - a) * x + 0.5 * (b + a)))
        return np.float128(0.5 * (b - a)*ans)

    def simp(self, a, b, f, n, split):

        dx = (b - a) / n
        X = np.linspace(a, b, split)
        Y = [f(x) for x in X]
        sum = dx / 3 * np.sum(Y[0:-1:2] + 4 * Y[1::2] + Y[2::2])
        return np.float32(sum)
        pass

    def integrate(self, f: callable, a: float, b: float, n: int) -> np.float32:
        """
        Integrate the function f in the closed range [a,b] using at most n 
        points. Your main objective is minimizing the integration error. 
        Your secondary objective is minimizing the running time. The assignment
        will be tested on variety of different functions. 
        
        Integration error will be measured compared to the actual value of the 
        definite integral. 
        
        Note: It is forbidden to call f more than n times. 
        
        Parameters
        ----------
        f : callable. it is the given function
        a : float
            beginning of the integration range.
        b : float
            end of the integration range.
        n : int
            maximal number of points to use.

        Returns
        -------
        np.float32
            The definite integral of f between a and b
        """

        if n == 2:
            return (f(b) - f(a)) / n
        if n % 2 != 0:
            split = n
            n = n - 1
        else:
            split = n - 1
            n = n - 1
        return np.float32(self.gauss(a, b, f, n))

        # replace this line with your solution

    def areabetween(self, f1: callable, f2: callable) -> np.float32:
        """
        Finds the area enclosed between two functions. This method finds 
        all intersection points between the two functions to work correctly. 
        
        Example: https://www.wolframalpha.com/input/?i=area+between+the+curves+y%3D1-2x%5E2%2Bx%5E3+and+y%3Dx

        Note, there is no such thing as negative area. 
        
        In order to find the enclosed area the given functions must intersect 
        in at least two points. If the functions do not intersect or intersect 
        in less than two points this function returns NaN.  
        This function may not work correctly if there is infinite number of 
        intersection points. 
        

        Parameters
        ----------
        f1,f2 : callable. These are the given functions

        Returns
        -------
        np.float32
            The area between function and the X axis

        """
        def creat_condition(f,a):
            if f(a) < 0:
                return lambda b: f(b) > 0
            else:
                return lambda b: f(b) < 0

        def falsePosition(f,x0, x1, e):
            x2 = x0
            while abs(f(x2)) >= e:
                x2 = x0 - (x1 - x0) * f(x0) / (f(x1) - f(x0))
                if f(x0) * f(x2) < 0:
                    x1 = x2
                else:
                    x0 = x2
            return x2

        def next_root(f,left, b, maxerr,to_add):
            right = left + to_add
            if (right > b):
                return None
            condition = creat_condition(f,left)
            div = 9
            while not condition(right):
                if abs(f(left)) < maxerr:
                    return left
                to_add = maxerr/div
                if div < 3:
                    div = random.randint(1,11)
                else:
                    div = div - 2
                left = left + to_add
                right = right + to_add
                if (right > b):
                    return None
            return falsePosition(f,left, right, maxerr)

        def inter_with_window(f,a, b, maxerr):
            ans = []
            eps = sys.float_info.epsilon
            x = next_root(f,a, b, maxerr/10,maxerr/10)
            while (x is not None):
                if x <= b:
                    ans.append(x)
                to_add = maxerr/10
                a = x + maxerr + eps
                x = next_root(f,a, b, maxerr,to_add)
            if len(ans) != 0 and abs(ans[len(ans) - 1] - b) < maxerr and abs(f(b)) < maxerr:
                ans.append(b)
            return ans
        # replace this line with your solution
        result = np.float128(0)
        f = lambda x: f1(x) - f2(x)
        xs = inter_with_window(f, 1, 100,0.01)
        n = 300
        length = len(xs)
        if length < 2:
            return np.float128(0)
        for i in range(length - 1):
            result = np.float128(result + np.abs(self.gauss(xs[i], xs[i + 1], f,  n)))
        return np.float32(result)


##########################################################################


import unittest
from sampleFunctions import *


class TestAssignment3(unittest.TestCase):

    def test_integrate_hard_case(self):
        ass3 = Assignment3()
        f1 = strong_oscilations()
        r = ass3.integrate(f1, 0.09, 10, 211)
        true_result = -7.78662 * 10 ** 33
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_integrate_float32(self):
        ass3 = Assignment3()
        f1 = RESTRICT_INVOCATIONS(11)(np.poly1d([-1, 0, 1]))
        r = ass3.integrate(f1, -1, 1, 10)
        self.assertEquals(r.dtype, np.float32)

    def test_polynomial_and_ln_case(self):
        ass3 = Assignment3()
        f = RESTRICT_INVOCATIONS(18)(lambda x: (x - 4) * (x - 2.5) - math.log(x, math.e))
        r = ass3.integrate(f, 1, 5.5, 18)
        true_result = - (88 * math.log((11 / 2), math.e) - 153) / 16
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_case_constant(self):
        ass3 = Assignment3()
        f = RESTRICT_INVOCATIONS(7)(lambda x: 5)
        r = ass3.integrate(f, -5, 10, 7)
        true_result = 75
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_polynomial(self):
        ass3 = Assignment3()
        f = RESTRICT_INVOCATIONS(8)(lambda x: x * x - 3 * x + 5)
        r = ass3.integrate(f, -0.07, 7.1, 8)
        true_result = 79.546131
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_sin1(self):
        ass3 = Assignment3()
        f = lambda x: math.sin((x * x))
        r = ass3.integrate(f, -4.4, 7.31, 40)
        true_result = 1.221308543955602
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_exp1(self):
        ass3 = Assignment3()
        f = lambda x: math.e ** (-2 * x * x)
        r = ass3.integrate(f, -2.71, 10.42, 19)
        true_result = 1.253314099967343
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_arctan(self):
        ass3 = Assignment3()
        f = lambda x: math.atan(x)
        r = ass3.integrate(f, -0.27, 5.63, 15)
        true_result = 6.074245208223067
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_sin2(self):
        ass3 = Assignment3()

        def f(x):
            return np.sin(x) / x

        r = ass3.integrate(f, -10.2, 8.4123, 20)
        true_result = 3.26705
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_ln1(self):
        ass3 = Assignment3()
        f = lambda x: 1 / math.log(x, math.e)
        r = ass3.integrate(f, 2.31, 9.74, 13)
        true_result = 4.60146
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_ln2(self):
        ass3 = Assignment3()
        f = lambda x: math.log(math.log(x, math.e), math.e)
        r = ass3.integrate(f, 2, 10, 15)
        true_result = 3.95291
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_sinln(self):
        ass3 = Assignment3()
        f = lambda x: math.sin(math.log(x, math.e))
        r = ass3.integrate(f, 3, 7, 13)
        true_result = 3.8853
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_strong_oscilations(self):
        ass3 = Assignment3()
        f = sampleFunctions.strong_oscilations()
        r = ass3.integrate(f, 1, 3, 9)
        true_result = 1.36924843371
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_area_between_3(self):
        ass3 = Assignment3()
        f1, f2 = lambda x: (x - 2) * (x - 5) * (x - 10), lambda x: (x - 2) * (x - 11)
        r = ass3.areabetween(f1, f2)
        true_result = 207 / 4 + 36 * math.sqrt(3)
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_area_between_4(self):
        ass3 = Assignment3()
        f1, f2 = lambda x: (x - 2) * (x - 5) * (x - 10) * (x - 15), lambda x: (x - 2) * (x - 11)
        r = ass3.areabetween(f1, f2)
        true_result = 2986.42
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_area_between_5(self):
        ass3 = Assignment3()
        f1, f2 = lambda x: (math.e ** x) / 37, lambda x: (x - 3) * (x - 4) * (x - 7)
        r = ass3.areabetween(f1, f2)
        true_result = 0.0095818
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_area_between_6(self):
        ass3 = Assignment3()
        f1, f2 = lambda x: (x - 2) * (x - 5) * (x - 10) * (x - 15), lambda x: (x - 2) * (x - 11)
        r = ass3.areabetween(f1, f2)
        true_result = 2986.42
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))

    def test_area_between_7(self):
        ass3 = Assignment3()
        f1, f2 = lambda x: (math.e ** x) / 4, lambda x: 5 * (x - 2) * (x - 3) * (x - 3.5) * (x - 6) * (x - 20)
        r = ass3.areabetween(f1, f2)
        true_result = 967.785
        self.assertGreaterEqual(0.001, abs((r - true_result) / true_result))


if __name__ == "__main__":
    unittest.main()
