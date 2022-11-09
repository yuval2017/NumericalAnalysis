"""
In this assignment you should fit a model function of your choice to data
that you sample from a given function.

The sampled data is very noisy so you should minimize the mean least squares
between the model you fit and the data points you sample.

During the testing of this assignment running time will be constrained. You
receive the maximal running time as an argument for the fitting method. You
must make sure that the fitting function returns at most 5 seconds after the
allowed running time elapses. If you take an iterative approach and know that
your iterations may take more than 1-2 seconds break out of any optimization
loops you have ahead of time.

Note: You are NOT allowed to use any numeric optimization libraries and tools
for solving this assignment.

"""
from functools import cmp_to_key
from math import hypot
import random
import scipy.optimize
import assignment1


class Assignment4A:
    def __init__(self):
        """
        Here goes any one time calculation that need to be made before
        solving the assignment for specific functions.
        """

        pass

    def get_cnrl_pointts(self, points, M_inv, T):
        x2 = self.inv_matrix(np.dot(np.transpose(T), T))
        x3 = np.dot(M_inv, x2)
        x4 = np.dot(x3, np.transpose(T))
        return np.dot(x4, points)

    def inv_matrix(self, A):
        n = len(A)
        AM = np.array(A, dtype=np.longdouble)
        I_M = np.identity(n, dtype=np.longdouble)
        ind = list(range(n))
        for k in range(n):
            fd_scaler = 1.0 / AM[k][k]
            for j in range(n):
                AM[k][j] *= fd_scaler
                I_M[k][j] *= fd_scaler
            for i in ind[0:k] + ind[k + 1:]:
                crScaler = AM[i][k]
                for j in range(n):
                    AM[i][j] = AM[i][j] - crScaler * AM[k][j]
                    I_M[i][j] = I_M[i][j] - crScaler * I_M[k][j]
        return I_M

    def take_sample(self, f, a, b):
        p = random.uniform(a, b)
        return p, f(p)

    def get_sorted_samples(self, f, a, b, n, to_left, maxtime, T_initial, ans):
        if (len(ans) != 0):
            points = ans
        else:
            points = [(a, f(a))]
        fb = f(b)
        for i in range(n - 2):
            if time.time() - T_initial >= maxtime - to_left:
                break
            points.append(self.take_sample(f, a, b))
        points.append((b, fb))
        return np.array(sorted(points, key=cmp_to_key(lambda p1, p2: p1[0] - p2[0])), dtype=np.longdouble)

    def get_distance_array(self, points):
        n = len(points)
        res = np.array([0] * n, dtype=np.longdouble)
        for i in range(1, n):
            point = points[i]
            prev = points[i - 1]
            res[i] = (res[i - 1] + (hypot(prev[0] - point[0], prev[1] - point[1])))
        return res

    def gen(self, div):
        def help_for_map2(val):
            val = val / div
            return [val ** 3, val ** 2, val ** 1, val ** 0]

        return help_for_map2

    def fit(self, f: callable, a: float, b: float, d: int, maxtime: float) -> callable:
        """
        Build a function that accurately fits the noisy data points sampled from
        some closed shape.

        Parameters
        ----------
        f : callable.
            A function which returns an approximate (noisy) Y value given X.
        a: float
            Start of the fitting range
        b: float
            End of the fitting range
        d: int
            The expected degree of a polynomial matching f
        maxtime : float
            This function returns after at most maxtime seconds.

        Returns
        -------
        a function:float->float that fits f between a and b
        """

        T_initial = time.time()
        if 10 <= maxtime < 15:
            points = self.get_sorted_samples(f, a, b, 400000, 4, maxtime, T_initial, [])
        elif 15 <= maxtime < 20:
            points = self.get_sorted_samples(f, a, b, 600000, 5, maxtime, T_initial, [])
        elif maxtime >= 20:
            points = self.get_sorted_samples(f, a, b, 1000000, 6, maxtime, T_initial, [])
        else:
            points = self.get_sorted_samples(f, a, b, 300000, 2, maxtime, T_initial, [])

        d_array = self.get_distance_array(points)

        div = d_array[len(d_array) - 1]
        func = self.gen(div)
        T = np.array([func(x) for x in d_array], dtype=np.longdouble)
        M_inv = np.array(
            [[0, 0, 0, 1],
             [0, 0, 1 / 3, 1],
             [0, 1 / 3, 2 / 3, 1],
             [1, 1, 1, 1]]
            , dtype=np.longdouble)
        """
        def starter(X, fcs):
            n = len(X) - 1
            def findVal(val):
                if val <= X[0]:
                    return fcs[0](val)[1]
                elif val > X[n - 1]:
                    return fcs[n - 2](val)[1]
                for j in range(0, n - 1):
                    if (val >= X[j]) and (val <= X[j + 1]):
                        t = (val - X[j]) / (X[j + 1] - X[j])
                        func = fcs[j]
                        ans = func(t)
                        res = ans[1]
                        return res

            return findVal
            """

        ass1 = assignment1.Assignment1()
        C = self.get_cnrl_pointts(points, M_inv, T)
        funcs = ass1.get_cubic(C[0], C[1], C[2], C[3])

        def callie_d(x):
            polier = ass1.get_coef(C[0][0], C[1][0], C[2][0], C[3][0])(x)
            polier = np.poly1d(polier)
            rts = scipy.optimize.fsolve(polier, np.array([0, 1]))
            real_r = 0
            for r in rts:
                if np.isreal(r):
                    real_r = r
                    break
            x, y = funcs(real_r)
            return y

        return lambda x: callie_d(x)


##########################################################################


from sampleFunctions import *


class TestAssignment4(unittest.TestCase):

    def test_return(self):
        f = NOISY(0.01)(poly(1, 1, 1))
        ass4 = Assignment4A()
        T = time.time()
        shape = ass4.fit(f=f, a=0, b=1, d=10, maxtime=5)
        T = time.time() - T
        self.assertLessEqual(T, 5)

    def test_delay(self):
        f = DELAYED(7)(NOISY(0.01)(poly(1, 1, 1)))

        ass4 = Assignment4A()
        T = time.time()
        shape = ass4.fit(f=f, a=0, b=1, d=10, maxtime=5)
        T = time.time() - T
        self.assertGreaterEqual(T, 5)

    def test_err(self):
        f = poly(1, 1, 1)
        nf = NOISY(1)(f)
        ass4 = Assignment4A()
        T = time.time()
        ff = ass4.fit(f=nf, a=0, b=1, d=10, maxtime=5)
        T = time.time() - T
        mse = 0
        for x in np.linspace(0, 1, 1000):
            self.assertNotEquals(f(x), nf(x))
            mse += (f(x) - ff(x)) ** 2
        mse = mse / 1000
        print(mse)


if __name__ == "__main__":
    unittest.main()
