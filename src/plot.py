"""Plotting functions"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import cmcrameri

import utils

DEFAULT_CMAP = cmcrameri.cm.batlow

units = {
    'u': '[$m.s^{-1}$]',
    'v': '[$m.s^{-1}$]',
    'h': '[$m$]',
    'f': '[$s^{-1}$]'
}


def fmt(x):
    """Format to latitude."""
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} °"


def show_var(model, var, cmap=DEFAULT_CMAP, title=None, show_lat=False):
    """Plot a variable of the current state of the simulation."""
    mat = getattr(model, var)[:model.ny, :model.nx]
    fig, ax = plt.subplots()
    pcol = ax.pcolormesh(model.x_list[:-1], model.y_list[:-1], mat, cmap=cmap)
    if show_lat:
        cs = ax.contour(model.x_list[:-1], model.y_list[:-1], utils.co(model.colat[:model.ny, :model.nx]), colors='black', linestyles=':', linewidths=1, levels=4)
        ax.clabel(cs, cs.levels, inline=True, fmt=fmt, fontsize=8)
    ax.set_title(var if title is None else title)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    plt.colorbar(pcol, label=units[var] if var in units.keys() else '')
    return fig


def show_var_anim(model, var, n_frames=None, cmap=DEFAULT_CMAP, title=None, show_lat=False):
    data = getattr(model, var + '_data')
    data = data.sel(time=data.time <= model.timestep * model.dt)
    v_max = np.max(data)
    v_min = np.min(data)
    n = data.shape[0]
    if n_frames is None:
        n_frames = n
    
    def update(frame):
        k = int(frame * n / n_frames)
        plt.clf()
        plt.pcolormesh(data.x, data.y, data[k], cmap=cmap, vmax=v_max, vmin=v_min)
        plt.colorbar(label=var)
        plt.title(f't = {data.time[k] / 86400 :.0f} d')
        plt.xlabel('X')
        plt.ylabel('Y')

    return FuncAnimation(plt.gcf(), update, frames=n_frames, interval=200)
