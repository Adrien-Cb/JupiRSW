"""Plotting functions"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import cmcrameri
import warnings

import utils
from config import FFMPEG_PATH

plt.style.use('ggplot')

if FFMPEG_PATH is not None and FFMPEG_PATH != '':
    plt.rcParams['animation.ffmpeg_path'] = FFMPEG_PATH

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


def show_var_anim(model, var, n_frames=None, cmap=DEFAULT_CMAP, title=None, polar_grid=True, 
                  save_as=None, filename='', save_dpi=200, fps=24):
    if var == 'pv':
        try:
            data = getattr(model.data, var)
        except AttributeError:
            data = model.data.compute_pv(model.h_0, model.f)
    else:
        data = getattr(model.data, var)
    data = data.sel(time=data.time <= model.timestep * model.dt)
    v_max = np.max(data)
    v_min = np.min(data)
    n = data.shape[0]
    if n_frames is None:
        n_frames = n

    # Create plot
    fig, ax = plt.subplots(1, 1)
    pcol = ax.pcolormesh(data.x, data.y, data[0], cmap=cmap, vmax=v_max, vmin=v_min)
    fig.colorbar(pcol, label=var if title is None else title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    # Add grid
    if polar_grid:
        lw = 0.5
        ls = ':'
        color = 'black'
        x_center = (model.x_list[-2] + model.x_list[0]) * .5
        y_center = (model.y_list[-2] + model.y_list[0]) * .5
        ax.axvline(x_center, linewidth=lw, linestyle=ls, color=color)
        ax.axhline(y_center, linewidth=lw, linestyle=ls, color=color)
        cs = ax.contour(model.x_list[:-1], model.y_list[:-1], utils.co(model.colat[:model.ny, :model.nx]), 
                        colors=color, linestyles=ls, linewidths=lw, levels=4)
        ax.clabel(cs, cs.levels, inline=True, fontsize=8)
        ax.axis('off')
    
    # Update plot
    def update(frame):
        k = int(frame * n / n_frames)
        pcol.set_array(data[k].values.ravel())
        ax.set_title(f't = {data.time[k] / 86400 :.0f} d')

    anim = FuncAnimation(fig, update, frames=n_frames, interval=200)

    # Save
    if filename == '':
        filename = utils.generate_output_name('anim')
    if save_as in ('mp4', 'MP4'):
        ffwriter = FFMpegWriter()
        try:
            anim.save(filename.strip('.mp4') + '.mp4', writer=ffwriter, dpi=save_dpi)
        except (PermissionError, FileNotFoundError):
            utils.warn('Could not save as MP4, check that ffmpeg executable is at the path specified in the configuration file.')
    elif save_as in ('gif', 'GIF'):
        anim.save(filename.strip('.gif') + '.gif', dpi=save_dpi, fps=fps)
    elif save_as is not None:
        utils.warn('Invalid save format, it should be \'mp4\' or \'gif\'.')

    return anim
