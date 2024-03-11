"""Core of the simulation"""

import numpy as np


def ddx(z, dx):
    """Compute x derivative of array z."""
    return (z[:, 1:] - z[:, :-1]) / dx


def ddy(z, dy):
    """Compute y derivative of array z."""
    return (z[1:, :] - z[:-1, :]) / dy


def avx(z):
    """Reduce x length of z by calculating the average of consecutive cells."""
    return (z[:, 1:] + z[:, :-1]) * 0.5


def avy(z):
    """Reduce y length of z by calculating the average of consecutive cells."""
    return (z[1:, :] + z[:-1, :]) * 0.5


def curl(u, v, dx, dy):
    """Compute curl of u and v."""
    omega = np.zeros((v.shape[0], u.shape[1]))
    omega[1:-1, :] = ddy(u, dy)
    omega[:, 1:-1] -= ddx(v, dx)
    return omega


def rhs(state, f, g, h_0, dx, dy):
    """Compute the RHS of the RSW equations.
    Input: state is a list of three arrays
    Output: a list of three arrays with the tendencies"""
    h_, u_, v_ = state

    h_x = np.zeros_like(u_)  # h with the same shape as u
    h_y = np.zeros_like(v_)  # h with the same shape as v
    h_x[:, 1:-1] = avx(h_0 + h_)
    h_y[1:-1, :] = avy(h_0 + h_)

    d_h = -(ddx(h_x * u_, dx) + ddy(h_y * v_, dy))

    # Vorticity at corner cell
    omega = f + curl(u_, v_, dx, dy)

    # Bernoulli function
    b = g * h_ + 0.5 * (avx(u_**2) + avy(v_**2))

    du = np.zeros_like(u_)
    dv = np.zeros_like(v_)

    # Centered discretization
    du[:, 1:-1] = -ddx(b, dx) + avy(omega[:, 1:-1]) * avy(avx(v_))
    dv[1:-1, :] = -ddy(b, dy) - avx(omega[1:-1, :]) * avx(avy(u_))

    return d_h, du, dv


def rk3(h, u, v, dt, f, g, h_0, dx, dy):
    state = h, u, v
    ds0 = rhs(state, f, g, h_0, dx, dy)
    s0 = [x + y * dt for x, y in zip(state, ds0)]

    ds1 = rhs(s0, f, g, h_0, dx, dy)
    s1 = [x + (y + z) * (dt / 4) for x, y, z in zip(state, ds0, ds1)]

    ds2 = rhs(s1, f, g, h_0, dx, dy)
    state = [
        x + (w + y + 4 * z) * (dt / 6) for x, w, y, z in zip(state, ds0, ds1, ds2)
    ]
    return state


def pv(u, v, h_tot, f, dx, dy):
    """Compute potential vorticity.
    Here h_tot is the total thickness of the layer, i.e. h + h_0."""
    return (curl(u, v, dx, dy) + f)[:-1, :-1] / h_tot