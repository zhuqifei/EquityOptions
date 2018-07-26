from solver.Solver import Solver

from numpy.core import empty, clip, zeros, exp, sqrt, ceil
from scipy import sparse
import scipy.sparse.linalg.dsolve as linsolve


class BSSolver(Solver):
    """Finite-difference solver for the Black-Scholes PDE in its most basic form.
           The problem to solve is given by:

            Function
              f = f(t,x) over the rectangle 0 <= t <= T, Bl <= x <= Bu.
            PDE
              rf = df/dt + rx df/dx + (1/2)(sigma x)^2 d^f/dx^2
            Boundary conditions
              given on the three sides of the rectangle
              t = T; x = Bl; x = Bu

           where r, sigma, T, Bl, Bu are given parameters.
        """

    def __init__(self, r, sigma, T, Bl, Bu, Fl, Fu, Fp, m, n, isAmericanOpt):

        super().__init__()

        """Initialize the finite-difference solver.

         Parameters:
          r     - interest rate
          sigma - volatility
          T     - maturity time
          Bl    - stock price on lower boundary
          Bu    - stock price on upper boundary
          Fl    - value of option on lower boundary
          Fu    - value of option on upper boundary
          Fp    - pay-off at maturity
          m     - number of time steps to take when discretizing PDE
          n     - number of points in x (stock price) domain
                  when discretizing PDE;  does not include the boundary points
        """

        self.r = r;
        self.sigma = sigma;
        self.T = T

        self.Bl = Bl;
        self.Bu = Bu
        self.Fl = Fl;
        self.Fu = Fu

        self.m = m;
        self.n = n

        # Step sizes
        self.dt = float(T) / m
        self.dx = float(Bu - Bl) / (n + 1)
        self.xs = Bl / self.dx

        # Grid that will eventually contain the finite-difference solution
        self.u = empty((m + 1, n))
        self.u[0, :] = Fp  # initial condition
        self.isAmerican = isAmericanOpt
        self.intrinsicValue = self.u[0, :]

    def build_sparse_explicit(self, s):
        """(internal) Set up the sparse matrix system for the explicit method."""

        A = sparse.lil_matrix((s.n, s.n))

        for j in range(0, s.n):
            xd = j + 1 + s.xs
            ssxx = (s.sigma * xd) ** 2

            A[j, j] = 1.0 - s.dt * (ssxx + s.r)

            if j > 0:
                A[j, j - 1] = 0.5 * s.dt * (ssxx - s.r * xd)
            if j < s.n - 1:
                A[j, j + 1] = 0.5 * s.dt * (ssxx + s.r * xd)

        s.A = A.tocsr()

    def build_sparse_implicit(self, s):
        """(internal) Set up the sparse matrix system for the implicit method."""

        C = sparse.lil_matrix((s.n, s.n))

        for j in range(0, s.n):
            xd = j + 1 + s.xs
            ssxx = (s.sigma * xd) ** 2

            C[j, j] = 1.0 + s.dt * (ssxx + s.r)

            if j > 0:
                C[j, j - 1] = 0.5 * s.dt * (-ssxx + s.r * xd)
            if j < s.n - 1:
                C[j, j + 1] = 0.5 * s.dt * (-ssxx - s.r * xd)

                # Store matrix with sparse LU decomposition already performed
        s.C = linsolve.splu(sparse.lil_matrix.tocsc(C))

        # Buffer to store right-hand side of the linear system Cu = v
        s.v = empty((n,))

    def build_sparse_crank_nicolson(self, s):
        """(internal) Set up the sparse matrices for the Crank-Nicolson method. """

        A = sparse.lil_matrix((s.n, s.n))
        C = sparse.lil_matrix((s.n, s.n))

        for j in range(0, s.n):
            xd = j + 1 + s.xs
            ssxx = (s.sigma * xd) ** 2

            A[j, j] = 1.0 - 0.5 * s.dt * (ssxx + s.r)
            C[j, j] = 1.0 + 0.5 * s.dt * (ssxx + s.r)

            if j > 0:
                A[j, j - 1] = 0.25 * s.dt * (+ssxx - s.r * xd)
                C[j, j - 1] = 0.25 * s.dt * (-ssxx + s.r * xd)
            if j < s.n - 1:
                A[j, j + 1] = 0.25 * s.dt * (+ssxx + s.r * xd)
                C[j, j + 1] = 0.25 * s.dt * (-ssxx - s.r * xd)

        s.A = A.tocsr()
        s.C = linsolve.splu(sparse.lil_matrix.tocsc(C))  # perform sparse LU decomposition

        # Buffer to store right-hand side of the linear system Cu = v
        s.v = empty((self.n,))

    def time_step_explicit(self, s, i):
        """(internal) Solve the PDE for one time step using the explicit method."""

        # Perform the next time step
        s.u[i + 1, :] = s.A * s.u[i, :]

        # and mix in the two other boundary conditions not accounted for above
        xdl = 1 + s.xs;
        xdu = s.n + s.xs
        s.u[i + 1, 0] += s.Fl[i] * 0.5 * s.dt * ((s.sigma * xdl) ** 2 - s.r * xdl)
        s.u[i + 1, s.n - 1] += s.Fu[i] * 0.5 * s.dt * ((s.sigma * xdu) ** 2 + s.r * xdu)

    def time_step_implicit(self, s, i):
        """(internal) Solve the PDE for one time step using the implicit method."""

        s.v[:] = s.u[i, :]

        # Add in the two other boundary conditions
        xdl = 1 + s.xs;
        xdu = s.n + s.xs
        s.v[0] -= s.Fl[i + 1] * 0.5 * s.dt * (-(s.sigma * xdl) ** 2 + s.r * xdl)
        s.v[s.n - 1] -= s.Fu[i + 1] * 0.5 * s.dt * (-(s.sigma * xdu) ** 2 - s.r * xdu)

        # Perform the next time step
        s.u[i + 1, :] = s.C.solve(s.v)

    def time_step_crank_nicolson(self, s, i):
        """(internal) Solve the PDE for one time step using the Crank-Nicolson method."""

        # Perform explicit part of time step
        s.v[:] = s.A * s.u[i, :]

        # Add in the two other boundary conditions
        xdl = 1 + s.xs;
        xdu = s.n + s.xs
        s.v[0] += s.Fl[i] * 0.25 * s.dt * (+(s.sigma * xdl) ** 2 - s.r * xdl)
        s.v[s.n - 1] += s.Fu[i] * 0.25 * s.dt * (+(s.sigma * xdu) ** 2 + s.r * xdu)
        s.v[0] -= s.Fl[i + 1] * 0.25 * s.dt * (-(s.sigma * xdl) ** 2 + s.r * xdl)
        s.v[s.n - 1] -= s.Fu[i + 1] * 0.25 * s.dt * (-(s.sigma * xdu) ** 2 - s.r * xdu)

        # Perform implicit part of time step
        s.u[i + 1, :] = s.C.solve(s.v)

    def solve(self, method='crank-nicolson'):
        """Solve the Black-Scholes PDE with the parameters given at initialization.

          Arguments:
           method - Indicates which finite-difference method to use.  Choices:
                    'explicit': explicit method
                    'implicit': implicit method
                    'crank-nicolson': Crank-Nicolson method
                    'smoothed-crank-nicolson':
                       Crank-Nicolson method with initial smoothing
                       by the implicit method
        """

        i_start = 0

        if method == 'implicit':
            self.build_sparse_implicit()
            time_step = self.time_step_implicit
        elif method == 'explicit':
            self.build_sparse_explicit()
            time_step = self.time_step_explicit
        elif method == 'crank-nicolson' or method is None:
            self.build_sparse_crank_nicolson()
            time_step = self.time_step_crank_nicolson
        elif method == 'smoothed-crank-nicolson':
            self.build_sparse_implicit()
            for i in range(0, 10):
                self.time_step_implicit(i)
            i_start = 10
            self.build_sparse_crank_nicolson()
            time_step = self.time_step_crank_nicolson
        else:
            raise ValueError('incorrect value for method argument')

        for i in range(i_start, m):
            time_step(i)

        return self.u
