#!/usr/bin/env python

# To run this, you need Python 2.x, and install the
# numpy and matplotlib packages for that language.

from pricer.Pricer import Pricer
from solver.BSSolver import BSSolver

from numpy.core import empty, clip, zeros, exp, sqrt, ceil
from numpy import linspace
from scipy import sparse
import scipy.sparse.linalg.dsolve as linsolve
import pylab
import matplotlib as matplot
import mpl_toolkits.mplot3d.axes3d as ax3d


class BSPricer(Pricer):

    def price(self):
        return NotImplemented

    def european_put(self, r, sigma, T, Bu, m, n, Bl=0.0, barrier=None, method=None, isAmerican=False):
        """Compute prices for a European-style put option."""

        X = linspace(0.0, B, n + 2)
        X = X[1:-1]

        Fp = clip(K - X, 0.0, 1e600)

        if barrier is None:
            Fu = zeros((m + 1,))
            Fl = K * exp(-r * linspace(0.0, T, m + 1))
        elif barrier == 'up-and-out':
            Fu = Fl = zeros((m + 1,))

        bss = BSSolver(r, sigma, T, Bl, Bu, Fl, Fu, Fp, m, n, isAmerican)
        return X, bss.solve(method)

    def european_call(self, r, sigma, T, Bu, m, n, Bl=0.0, barrier=None, method=None, isAmerican=False):
        """Compute prices for a European-style call option."""

        X = linspace(0.0, B, n + 2)
        X = X[1:-1]

        Fp = clip(X - K, 0.0, 1e600)

        if barrier is None:
            Fu = B - K * exp(-r * linspace(0.0, T, m + 1))
            Fl = zeros((m + 1,))
        elif barrier == 'up-and-out':
            Fu = Fl = zeros((m + 1,))

        bss = BSSolver(r, sigma, T, Bl, Bu, Fl, Fu, Fp, m, n, isAmerican)
        return X, bss.solve(method)

    def plot_solution(self, T, X, u):
        # The surface plot
        Xm, Tm = pylab.meshgrid(X, linspace(T, 0.0, u.shape[0]))
        fig_surface = pylab.figure()
        ax = ax3d.Axes3D(fig_surface)
        ax.plot_surface(Xm, Tm, u)
        ax.set_ylabel('Time $t$')
        ax.set_xlabel('Stock price $x$')
        ax.set_zlabel('Option value $f(t,x)$')

        # The color temperature plot
        fig_color = pylab.figure()
        ax = pylab.gca()
        ax.set_xlabel('Time $t$')
        ax.set_ylabel('Stock price $x$')
        ax.imshow(u.T, interpolation='bilinear', origin='lower',
                  cmap=matplot.cm.hot, aspect='auto', extent=(T, 0.0, X[0], X[-1]))

        # Plot of price function at time 0
        fig_zero = pylab.figure()
        pylab.plot(X, u[m - 1, :])
        ax = pylab.gca()
        ax.set_xlabel('Stock price $x$')
        ax.set_ylabel('Option value $f(0,x)$')

        return fig_surface, fig_color, fig_zero


if __name__ == "__main__":
    def parse_options():
        from optparse import OptionParser

        parser = OptionParser()

        parser.add_option("-r", "--interest", dest="r", default="0.10",
                          help="interest rate [default: %default]")
        parser.add_option("-v", "--volatility", dest="sigma", default="0.40",
                          help="volatility [default: %default]")

        parser.add_option("-K", "--strike", dest="K", default="50.00",
                          help="strike price [default: %default]")
        parser.add_option("-T", "--maturity", dest="T", default="0.5",
                          help="maturity time [default: %default]")
        parser.add_option("-B", "--bound", dest="B", default="100.00",
                          help="upper bound on stock price [default: %default]")

        parser.add_option("-m", "--time-steps", type="int", dest="m", default="100",
                          help="number of time steps [default: %default]")
        parser.add_option("-n", "--space-steps", dest="n", default="200",
                          help="number of steps in stock-price space [default: %default]")
        parser.add_option("--dt", dest="dt", type="float", help="time step size")
        parser.add_option("--dx", dest="dx", type="float", help="stock-price step size")
        parser.add_option("--method", dest="method", help="finite-difference method")

        parser.add_option("-C", "--call", dest="is_call", action="store_true",
                          help="value a European-style call option")
        parser.add_option("-P", "--put", dest="is_put", action="store_true",
                          help="value a European-style put option")
        parser.add_option("--barrier", dest="barrier",
                          help="value a barrier option")

        parser.add_option("--plot", dest="plot", action="store_true",
                          help="plot results")
        parser.add_option("--save-plot", dest="save_plot", action="store_true",
                          help="save plots to EPS files")

        (options, args) = parser.parse_args()
        return options


    options = parse_options()

    # Parameters
    r = float(options.r)
    sigma = float(options.sigma)
    K = float(options.K)
    T = float(options.T)
    B = float(options.B)
    dtt = float(options.dt)
    isAmericanOption = True

    m = int(options.m)
    n = int(options.n)

    if options.dt is not None:
        m = int(ceil(T / float(options.dt)))
    if options.dx is not None:
        n = int(ceil(B / float(options.dx))) - 1

    if options.is_put:
        X, u = BSPricer.european_put(r, sigma, T, B, m, n,
                                     barrier=options.barrier, method=options.method, isAmerican=isAmericanOption)
    else:
        X, u = BSPricer().european_call(r, sigma, T, B, m, n,
                             barrier=options.barrier, method=options.method, isAmerican=isAmericanOption)

    # Print out results at time 0
    print("Stock price x    Option price f(0,x)")
    for i in range(0, n):
        print("%10.4f         %10.4f " % (X[i], u[m, i]))

    # Generate plots if user requests
    '''if options.plot or options.save_plot:

        golden_mean = (sqrt(5.0) - 1.0) / 2.0
        pylab.rcParams.update( \
            {'backend': 'ps',
             'ps.usedistiller': 'xpdf',
             'axes.labelsize': 10,
             'text.fontsize': 10,
             'xtick.labelsize': 8,
             'ytick.labelsize': 8,
             'figure.figsize': [7.0, golden_mean * 7.0],
             'text.usetex': True})

        fig_surface, fig_color, fig_zero = plot_solution(T, X, u)

        # Show figures
        #if options.plot:
        #    pylab.show()

        # Save figures to EPS formats
        if options.save_plot:
            fig_surface.savefig('bs-surface.eps')
            #fig_color.savefig('bs-color.eps')
            fig_zero.savefig('bs-zero.eps')
        '''
