"""Rotating shallow water model centered on the pole"""

import numpy as np
import xarray as xr
from scipy.special import gammaincc, gamma

from config import *
from data import Data
import utils
import core


def h_vort(r, h_0, r_m, ro, bu, b):
    """Return a vortex h at a distance r from the center"""
    x_ = 1/b * (r/r_m)**b
    s_ = 2/b
    gam = gamma(s_) * gammaincc(s_, x_)
    return + h_0 * (ro/bu * np.exp(1/b) * np.exp(2/b - 1) * gam)


def v_vort(r, v_m, r_m, b):
    """Return a vortex azimuthal velocity at a distance r from the center."""
    return - v_m * r/r_m * np.exp(1/b * (1 - (r/r_m)**b))


def ini_point(x, y, vort_centers, r_m, v_m, b, h_0, ro, bu):
    '''Return initial u, v, h at given point'''
    h = 0
    u = 0
    v = 0
    for x0, y0 in vort_centers:
        r_ = utils.r(x - x0, y - y0)
        h += h_vort(r_, h_0, r_m, ro, bu, b)
        if r_ < 1e-12:
            continue
        v_az = v_vort(r_, v_m, r_m, b)
        u -= v_az * (y-y0) / r_
        v += v_az * (x-x0) / r_
    return u, v, h


def make_vort_centers(vort_lat, vort_number, r_max, lat_min):
    """Return vortex centers.
    `vort_lat` and `vort_number` are either two scalar or two iterables of same length.
    """
    try:
        vort_lat[0]
    except TypeError:
        vort_lat = [vort_lat]
        vort_number = [vort_number]
    centers = []
    for lat, number in zip(vort_lat, vort_number):
        r_ = utils.co(lat) / utils.co(lat_min) * r_max
        angle_step = 2 * np.pi / number 
        for k in range(number):
            centers.append((r_ * np.cos(k*angle_step), r_ * np.sin(k*angle_step)))
    return centers


class Model():
    """Rotating shallow water model"""
    def __init__(self, nx, ny=None, 
                 lat_min=LAT_MIN, lat_sponge=LAT_SPONGE,
                 r_planet=R, t_planet=T, g=G,
                 ro=RO, bu=BU, b=B, r_m=R_M, cfl=CFL, output_dt=OUTPUT_TIMESTEP):
        self.nx = nx
        self.ny = ny if ny is not None else nx
        self.timestep = 0

        self.lat_min = lat_min
        self.lat_sponge = lat_sponge

        self.r_max = utils.co(LAT_MIN) / 180 * np.pi * r_planet
        self.dx = 2*self.r_max / (self.nx-1)
        self.dy = 2*self.r_max / (self.ny-1)

        self.bu = bu
        self.ro = ro
        self.b = b
        self.r_m = r_m
        self.g = g

        self.f_0 = 4*np.pi/t_planet
        self.v_m = self.ro * self.f_0 * self.r_m
        self.h_0 = self.bu * (self.f_0 * self.r_m)**2 / self.g

        u_max = np.sqrt(self.g * self.h_0)
        self.dt = int(cfl / np.sqrt((u_max/self.dx)**2 + (u_max/self.dy)**2))

        output_dt = utils.parse_time(output_dt)
        self.output_nt = max(int(output_dt / self.dt), 1)

        # Coordinates
        self.x_list = np.linspace(-self.r_max, self.r_max + self.dx, self.nx + 1)
        self.y_list = np.linspace(-self.r_max, self.r_max + self.dy, self.ny + 1)
        self.x, self.y = np.meshgrid(self.x_list, self.y_list)
        self.r = utils.r(self.x, self.y)
        self.colat = utils.colat(self.x, self.y, self.lat_min, self.r_max)

        # Prognostic variables
        self.u = np.zeros((self.ny, self.nx + 1))
        self.v = np.zeros((self.ny + 1, self.nx))
        self.h = np.zeros((self.ny, self.nx))

        # Coriolis parameter
        self.f = self.f_0 * np.cos(np.pi * self.colat / 180)

        # Sponge coefficient
        self.sponge = np.maximum((utils.co(self.lat_min) - np.maximum(self.colat, utils.co(self.lat_sponge)))
                                 / (utils.co(self.lat_min) - utils.co(self.lat_sponge)), 0)
        
        # Storage
        self.data = Data(self.x_list, self.y_list, self.output_nt * self.dt)
        
        utils.log(f'Configuration successfully created. Timestep: dt = {self.dt} s')

    def initialize(self, vort_lat, vort_number):
        """Set initial conditions.
        `vort_lat` and `vort_number` are either two scalars or two iterables of same length.
        """

        # Get list of vortex centers
        centers = make_vort_centers(vort_lat, vort_number, self.r_max, self.lat_min)

        # Initialize u, v and h
        for i in range(self.nx):
            for j in range(self.ny):
                self.u[j, i], self.v[j, i], self.h[j, i] = ini_point(
                    self.x[j, i], self.y[j, i], 
                    vort_centers=centers,
                    r_m=self.r_m, v_m=self.v_m, b=self.b, h_0=self.h_0, ro=self.ro, bu=self.bu 
                )

        # Initialize dataset (for storing the evolution of the system)
        self.data.store_state(0, self.u, self.v, self.h)

    def step(self):
        """Run one step of the simulation"""
        h, u, v = core.rk3(self.h, self.u, self.v,
                           dt=self.dt, f=self.f, g=self.g, h_0=self.h_0,
                           dx=self.dx, dy=self.dy)
        self.h = h
        self.u = u * self.sponge[:-1, :]
        self.v = v * self.sponge[:, :-1]
        self.timestep += 1
        if self.timestep % self.output_nt == 0:
            self.data.store_state(self.timestep * self.dt, self.u, self.v, self.h)

    def run(self, time):
        """Run simulation for given time.
        `time` can be a duration in seconds or a string such as `'2d 1h 15m 45s'`."""
        time = utils.parse_time(time)
        t_0 = self.timestep * self.dt
        step_init = self.timestep
        while self.timestep * self.dt - t_0 < time:
            self.step()
        utils.log(f'Successfully ran the simulation for {self.timestep - step_init} timesteps ({utils.sec_to_str(time)}).')

    def save_nc(self, filename=''):
        """Save output as NETCDF file."""
        self.data.save_nc(filename=filename)
