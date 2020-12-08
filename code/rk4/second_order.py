"""
Module for solving second order diferential equations using Runge Kutta 4 method.
"""
import sys

import numpy as np


def rk4_step(ode, x, y, ydot, h, *args, **kwargs):
    """
    @returns y(x+h), ydot(x+h)
    """
    k1 = h * ydot
    l1 = h * ode(x, y, ydot, *args, **kwargs)

    k2 = h * (ydot + 0.5*l1)
    l2 = h * ode(x + 0.5*h, y + 0.5*k1, ydot + 0.5*l1, *args, **kwargs)

    k3 = h * (ydot + 0.5*l2)
    l3 = h * ode(x + 0.5*h, y + 0.5*k2, ydot + 0.5*l2, *args, **kwargs)

    k4 = h * (ydot + l3)
    l4 = h * ode(x + h, y + k3, ydot + l3, *args, **kwargs)

    return y + (k1 + 2*k2 + 2*k3 + k4)/6, ydot + (l1 + 2*l2 + 2*l3 + l4)/6


def solver(ode, y0, ydot0, xmin, xmax, dx, *args, term=None, **kwargs):
    """
    @param ode: function that returns an array with all the derivatives. Must be solution to
        d2y/dx2 = ode(x, y, ydot)
    @param y0: an array of the same size with the initial values.
    @param ydot0: an array of the same size with the initial values for the derivatives.
    @param xmin: ode's domain minimum
    @param xmax: ode's domain maximum
    @param dx: step size
    @param term: function that evaluates termination condition

    @params *args, *kwargs: passed to ode

    side note: x is the independent variable for all.
    """
    if bool(term):
        if not callable(term):
            raise TypeError("term must be a function")
    if len(y0) != len(ydot0):
        raise ValueError("y0 and ydot0 must have same size")

    N = len(y0)
    nsteps = int((xmax - xmin)/dx)

    x = np.linspace(xmin, xmax, nsteps)
    u = np.zeros([nsteps, N])
    y = np.zeros([nsteps, N])

    y[0] = y0
    u[0] = ydot0

    for i in range(nsteps - 1):
        y[i+1], u[i+1] = rk4_step(ode, x[i], y[i], u[i], dx, *args, **kwargs)

        if bool(term):
            if term(y[i+1]):
                nsteps = i + 1
                break

    return x[:nsteps], y[:nsteps, :]


def n_solver(ode, tmin, tmax, x0, v0, N):
    """
    Solves a second order diff equation for a system of n bodies

    ode must return x_i(t+dt), v_i(t+dt), where i is the body id and x & v are three dimensional arrays
    x0 initial positions for each body
    v0 initial velocities for each body
    N number of time steps

    Returns
    x: positions, where
        x[body id][time step][coordinate(x, y or z)]
    v: velocities, where
        v[body id][time step][coordinate(x, y or z)]
    """
    dt = (tmax-tmin)/N

    P = len(x0)  # number of bodies

    x = np.zeros((P, N, 3))  # body, time, axis
    v = np.zeros((P, N, 3))

    x[:, 0] = x0
    v[:, 0] = v0

    unitColor = '\033[5;36m\033[5;47m'
    endColor = '\033[0;0m\033[0;0m'
    count = N-1

    for t in range(count):
        incre = int(50.0 / count * t)
        sys.stdout.write('\r')
        if t != (count - 1):
            sys.stdout.write('|%s%s%s%s| %d%%' % (
                unitColor, '\033[7m' + ' '*incre + ' \033[27m', endColor, ' '*(49-incre), 2*incre))

        else:
            sys.stdout.write('|%s%s%s| %d%%' % (
                unitColor, '\033[7m' + ' '*20 + 'COMPLETE!' + ' '*21 + ' \033[27m', endColor, 100))

        sys.stdout.flush()
        # print(t, end='\r', flush=True) # I included a "timer", this can be computationally demanding and I like knowing the progress
        for i in range(P):
            x[i][t+1], v[i][t +
                            1] = np.array(rk4_step(ode, 0., x[:, t], v[:, t], dt, i, P))[:, i]

    return x, v
