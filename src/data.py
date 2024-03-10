"""Handle data"""

import xarray as xr
import numpy as np

import utils



class Data():
    """Store u, v and h over time"""
    def __init__(self, x_list, y_list, output_dt):
        self.x_list = x_list
        self.y_list = y_list
        self.output_dt = output_dt

        self.nx = x_list.size - 1
        self.ny = y_list.size - 1
        self.time = 0

        self.make_dataset()

    def make_dataset(self):
        """Create dataset of u, v, h"""
        u = xr.DataArray(np.full((1, self.ny, self.nx + 1), np.nan), dims=('time', 'y', 'x'))
        v = xr.DataArray(np.full((1, self.ny + 1, self.nx), np.nan), dims=('time', 'y', 'x'))
        h = xr.DataArray(np.full((1, self.ny, self.nx), np.nan), dims=('time', 'y', 'x'))
        self.u = u.assign_coords(x=self.x_list, y=self.y_list[:-1], time=[0])
        self.v = v.assign_coords(x=self.x_list[:-1], y=self.y_list, time=[0])
        self.h = h.assign_coords(x=self.x_list[:-1], y=self.y_list[:-1], time=[0])

    def extend_dataset(self, n=100):
        """Extend dataset allocation (double the size or add n if size < n).
        This allows to indefinitely extend storage and thus simulation time while keeping a negligible
        amount of array extensions (which is costly due to reallocation)."""
        t_max = self.u.time.max()
        n = max(n, self.u.time.size)
        t_list = np.arange(t_max + self.output_dt, t_max + (n+1)*self.output_dt, self.output_dt)
        new_u = np.full((n, self.ny, self.nx + 1), np.nan)
        new_v = np.full((n, self.ny + 1, self.nx), np.nan)
        new_h = np.full((n, self.ny, self.nx), np.nan)
        new_u = xr.DataArray(new_u, dims=('time', 'y', 'x')).assign_coords(x=self.x_list, y=self.y_list[:-1], time=t_list)
        new_v = xr.DataArray(new_v, dims=('time', 'y', 'x')).assign_coords(x=self.x_list[:-1], y=self.y_list, time=t_list)
        new_h = xr.DataArray(new_h, dims=('time', 'y', 'x')).assign_coords(x=self.x_list[:-1], y=self.y_list[:-1], time=t_list)
        self.u = xr.concat([self.u, new_u], dim='time')
        self.v = xr.concat([self.v, new_v], dim='time')
        self.h = xr.concat([self.h, new_h], dim='time')
        utils.log(f'Dataset extension by {n} (total size of {self.u.time.size} timesteps).')
    
    def store_state(self, time, u, v, h):
        """Store current state of simulation"""
        try:
            d = dict(time=time)
            self.u.loc[d] = u
            self.v.loc[d] = v
            self.h.loc[d] = h
            self.time = max(self.time, time)
        except KeyError:
            self.extend_dataset()
            self.store_state(time, u, v, h)

    def save_nc(self, filename=''):
        """Save output as NETCDF file."""
        ds = xr.Dataset({'u': self.u, 'v': self.v, 'h': self.h})
        ds = ds.sel(time=ds.time <= self.time)
        if filename == '':
            filename = utils.generate_output_name('output')
        ds.to_netcdf(f'{filename}.nc')