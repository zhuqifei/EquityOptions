from solver.Solver import Solver

from scipy.optimize import optimize


class BSSolver(Solver):

    def solve(self, f, target):
        def f2(x):
            return abs(f(x) - target)
        return optimize.brent(f2)
