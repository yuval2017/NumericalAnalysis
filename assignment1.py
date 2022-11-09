"""
In this assignment you should interpolate the given function.
"""
import operator
from functools import reduce
import scipy.optimize


class Assignment1:
    def __init__(self):
        """
        Here goes any one time calculation that need to be made before 
        starting to interpolate arbitrary functions.
        """

        pass

    def do_chev(self):
        counter = 0

        def f():
            nonlocal counter
            if counter == 0:
                counter = counter + 1
                return True
            return False

        return f

    def interpolate_l(self, x_array, y_array):
        def g(x):
            def _basis(j):
                p = [(x - x_array[m]) / (x_array[j] - x_array[m]) for m in range(k) if m != j]
                return reduce(operator.mul, p)

            assert len(x_array) != 0 and (
                    len(x_array) == len(y_array))
            k = len(x_array)
            return sum(_basis(j) * y_array[j] for j in range(k))

        return g

    def TDMAsolver(self, a, b, c, d):
        nf = len(d)
        a_temp, b_temp, c_temp, d_temp = map(lambda arr: np.array(arr, dtype=np.longdouble), (a, b, c, d))
        for i in range(1, nf):
            m_c = (a_temp[i - 1] / b_temp[i - 1])
            b_temp[i] = (b_temp[i] - m_c * c_temp[i - 1])
            d_temp[i] = (d_temp[i] - m_c * d_temp[i - 1])

        ans = b_temp
        ans[-1] = (d_temp[-1] / b_temp[-1])
        for i in range(nf - 2, -1, -1):
            ans[i] = (d_temp[i] - c_temp[i] * ans[i + 1]) / b_temp[i]
        return ans

    def get_bezier_coef(self, points):
        n = len(points) - 1

        a_down = [1] * (n - 1)
        a_down[n - 2] = 2
        b_middle = [4] * (n)
        b_middle[n - 1] = 7
        b_middle[0] = 2
        c_up = [1] * (n - 1)

        P = [2 * (2 * points[i] + points[i + 1]) for i in range(n)]
        P[0] = points[0] + 2 * points[1]
        P[n - 1] = 8 * points[n - 1] + points[n]

        X, Y = [], []
        for p in P:
            X.append(p[0])
            Y.append(p[1])

        ptsX = self.TDMAsolver(a_down, b_middle, c_up, X)
        ptsY = self.TDMAsolver(a_down, b_middle, c_up, Y)
        ptsA = np.array([(x, y) for x, y in zip(ptsX, ptsY)], dtype=np.longdouble)

        B = [0] * n
        for i in range(n - 1):
            B[i] = 2 * points[i + 1] - ptsA[i + 1]
        B[n - 1] = (ptsA[n - 1] + points[n]) / 2

        return ptsA, B

    def get_coef(self, a, b, c, d):
        return lambda val: [(-a + 3 * b - 3 * c + d), (3 * a - 6 * b + 3 * c), (
                -3 * a + 3 * b), a - val]

    def get_cubic(self, a, b, c, d):
        return lambda t: np.array(
            np.power(1 - t, 3) * a + 3 * np.power(1 - t, 2) * t * b + 3 * (1 - t) * np.power(t, 2) * c + np.power(
                t, 3) * d, dtype=np.longdouble)

    def get_bezier_cubic(self, points):

        A, B = self.get_bezier_coef(points)
        ans1 = [
            self.get_cubic(points[i], A[i], B[i], points[i + 1])
            for i in range(len(points) - 1)
        ]
        ans2 = [self.get_coef(points[i][0], A[i][0], B[i][0], points[i + 1][0]) for i in range(len(points) - 1)]
        return ans1, ans2
    def find_root(self,pol, func):
        polier = np.poly1d(pol)
        rts = scipy.optimize.fsolve(polier, np.array([0, 1]))
        real_r = 0
        for r in rts:
            if np.isreal(r):
                real_r = r
        x, y = func(real_r)
        return y

    def findVal(self,val,funcs,xs,n,poliers, with_root = False):
        if val <= xs[0]:
            return funcs[0](val)[1]
        elif val > xs[n - 1]:
            return funcs[n - 2](val)[1]
        for j in range(0, n - 1):
            if (val >= xs[j]) and (val <= xs[j + 1]):
                if with_root:
                    yi = self.find_root(poliers[j](val), funcs[j])
                    return yi
                t = (val - xs[j]) / (xs[j + 1] - xs[j])
                func = funcs[j]
                xi, yi = func(t)
                if abs(xi - val) > 0.01:
                    yi = self.find_root(poliers[j](val), funcs[j])
                return yi

    def interpolate(self, f: callable, a: float, b: float, n: int) -> callable:
        """
        Interpolate the function f in the closed range [a,b] using at most n 
        points. Your main objective is minimizing the interpolation error.
        Your secondary objective is minimizing the running time. 
        The assignment will be tested on variety of different functions with 
        large n values. 
        
        Interpolation error will be measured as the average absolute error at 
        2*n random points between a and b. See test_with_poly() below. 

        Note: It is forbidden to call f more than n times. 

        Note: This assignment can be solved trivially with running time O(n^2)
        or it can be solved with running time of O(n) with some preprocessing.
        **Accurate O(n) solutions will receive higher grades.** 
        
        Note: sometimes you can get very accurate solutions with only few points, 
        significantly less than n. 
        
        Parameters
        ----------
        f : callable. it is the given function
        a : float
            beginning of the interpolation range.
        b : float
            end of the interpolation range.
        n : int
            maximal number of points to use.

        Returns
        -------
        The interpolating function.
        """

        # replace this line with your solution to pass the second test
        xs = np.linspace(a, b, n)
        points = np.array([(x, f(x)) for x in xs], dtype=np.longdouble)
        funcs, poliers = self.get_bezier_cubic(points)
        return lambda x: self.findVal(x, funcs, xs, n, poliers)
        """
        def cc():
            to_lag = {'value': 0}

            def chose_method(x):
                nonlocal to_lag
                to_lag['value'] += 1
                if to_lag['value'] > 200000:
                    return self.findVal(x,funcs,xs,n,poliers)
                else:
                    return lag(x)
            return chose_method
            """




##########################################################################

import unittest
from functionUtils import *
from tqdm import tqdm


class TestAssignment1(unittest.TestCase):
    def test_with_poly(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0

        d = 30
        num_of_points = 100
        for i in tqdm(range(100)):
            a = np.random.randn(d)

            f = np.poly1d(a)

            ff = ass1.interpolate(f, -10, 10, num_of_points)

            xs = np.random.random(2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("x^30 polynomial: " + str(T) + "[s]")
        print("x^30 polynomial: " + str(mean_err) + "[mean_err]")

    def test_with_poly_restrict(self):
        ass1 = Assignment1()
        a = np.random.randn(5)
        f = RESTRICT_INVOCATIONS(10)(np.poly1d(a))
        ff = ass1.interpolate(f, -10, 10, 10)
        xs = np.random.random(20)
        for x in xs:
            yy = ff(x)

    def test_with_sin(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 100
        for i in tqdm(range(100)):

            f = np.sin

            ff = ass1.interpolate(f, -10, 10, num_of_points)

            xs = np.random.uniform(low=-10, high=10, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("sin(x): " + str(T) + "[s]")
        print("sin(x): " + str(mean_err) + "[mean_err]")

    def test_with_y_5(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 100
        for i in tqdm(range(100)):

            f = lambda x: 5

            ff = ass1.interpolate(f, -10, 10, num_of_points)

            xs = np.random.uniform(low=-10, high=10, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("y=5: " + str(T) + "[s]")
        print("y=5: " + str(mean_err) + "[mean_err]")

    def test_with_sin_x_2(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 50
        for i in tqdm(range(100)):

            f = lambda x: np.sin(x ** 2)

            ff = ass1.interpolate(f, -1, 5, num_of_points)

            xs = np.random.uniform(low=-1, high=5, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("sin(x^2): " + str(T) + "[s]")
        print("sin(x^2): " + str(mean_err) + "[mean_err]")

    def test_with_e_with_exponent(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 100
        for i in tqdm(range(100)):

            f = lambda x: np.exp(-2 * (x ** 2))

            ff = ass1.interpolate(f, -2, 4, num_of_points)

            xs = np.random.uniform(low=-2, high=4, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("e^(-2x^2): " + str(T) + "[s]")
        print("e^(-2x^2): " + str(mean_err) + "[mean_err]")

    def test_with_arctan(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 100
        for i in tqdm(range(100)):

            f = lambda x: np.arctan(x)

            ff = ass1.interpolate(f, -5, 5, num_of_points)

            xs = np.random.uniform(low=-5, high=5, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("arctan: " + str(T) + "[s]")
        print("arctan: " + str(mean_err) + "[mean_err]")

    def test_with_sinx_div_x(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 100
        for i in tqdm(range(100)):

            f = lambda x: np.sin(x) / x

            ff = ass1.interpolate(f, 0.00001, 10, num_of_points)

            xs = np.random.uniform(low=0.00001, high=10, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("sin(x)/x: " + str(T) + "[s]")
        print("sin(x)/x: " + str(mean_err) + "[mean_err]")

    def test_with_1_div_lnx(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 200
        for i in tqdm(range(100)):

            f = lambda x: 1 / np.log(x)

            ff = ass1.interpolate(f, 0.00001, 0.9999, num_of_points)

            xs = np.random.uniform(low=0.00001, high=0.9999, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("1/ln(x): " + str(T) + "[s]")
        print("1/ln(x): " + str(mean_err) + "[mean_err]")

    def test_with_1_div_log(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 200
        for i in tqdm(range(100)):

            f = lambda x: pow(e, pow(e, x))

            ff = ass1.interpolate(f, 0.00001, 0.9999, num_of_points)

            xs = np.random.uniform(low=0.00001, high=0.9999, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("f9 " + str(T) + "[s]")
        print("f9: " + str(mean_err) + "[mean_err]")

    def test_with_e_e_x(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 200
        for i in tqdm(range(100)):

            f = lambda x: np.exp(np.exp(x))

            ff = ass1.interpolate(f, -2, 2, num_of_points)

            xs = np.random.uniform(low=-2, high=2, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("e^e^x: " + str(T) + "[s]")
        print("e^e^x: " + str(mean_err) + "[mean_err]")

    def test_with_ln_ln_x(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 100
        for i in tqdm(range(100)):

            f = lambda x: np.log(np.log(x))

            ff = ass1.interpolate(f, 1.00001, 30, num_of_points)

            xs = np.random.uniform(low=1.00001, high=30, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("ln(ln(x)): " + str(T) + "[s]")
        print("ln(ln(x)): " + str(mean_err) + "[mean_err]")

    def test_with_a_polynomial(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 100
        for i in tqdm(range(100)):

            f = lambda x: 5 * (x ** 2) - 10 * x + 1

            ff = ass1.interpolate(f, 3, 10, num_of_points)

            xs = np.random.uniform(low=3, high=10, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("5x^2-10x+1:" + str(T) + "[s]")
        print("5x^2-10x+1: " + str(mean_err) + "[mean_err]")

    def test_with_a_exp_sin(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 100
        for i in tqdm(range(100)):

            f = lambda x: 2 * (1 / (x * 2)) * np.sin(1 / x)

            ff = ass1.interpolate(f, 3, 10, num_of_points)

            xs = np.random.uniform(low=3, high=10, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("2^(1/x^2 )*sin(1/x):" + str(T) + "[s]")
        print("2^(1/x^2 )*sin(1/x): " + str(mean_err) + "[mean_err]")

    def test_with_e_ln(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 100
        for i in tqdm(range(100)):

            f = lambda x: np.log(np.exp(x))

            ff = ass1.interpolate(f, -10, 10, num_of_points)

            xs = np.random.uniform(low=-10, high=10, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("ln(e^x):" + str(T) + "[s]")
        print("ln(e^x): " + str(mean_err) + "[mean_err]")

    def test_with_ln_e(self):
        T = time.time()

        ass1 = Assignment1()
        mean_err = 0
        num_of_points = 100
        for i in tqdm(range(100)):

            f = lambda x: np.exp(np.log(x))

            ff = ass1.interpolate(f, 0.0000001, 10, num_of_points)

            xs = np.random.uniform(low=0.0000001, high=10, size=2 * num_of_points)
            err = 0
            for x in xs:
                yy = ff(x)
                y = f(x)
                err += abs(y - yy)

            err = err / (2 * num_of_points)
            mean_err += err
        mean_err = mean_err / 100

        T = time.time() - T
        print("e^(lnx):" + str(T) + "[s]")
        print("e^(lnx): " + str(mean_err) + "[mean_err]")


if __name__ == "main":
    unittest.main()
