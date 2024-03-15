"""Utility functions"""

import numpy as np
from datetime import datetime, timedelta
import re
from pathlib import Path

from config import LOG

def r(x, y):
    """Calculate distance of (x, y) to (0, 0)."""
    return np.sqrt(x**2 + y**2)


def co(lat):
    """Calculate colatitude from latitude or latitude from colatitude in degrees."""
    return 90 - lat


def colat(x, y, lat_min, r_max):
    '''Return colatitude of given x, y point'''
    return r(x, y) * co(lat_min) / r_max


def log(text):
    if not LOG:
        return
    print(f'[{datetime.now().strftime(r"%H:%M:%S")}] {text}')


def warn(text):
    log('WARNING | ' + text)


def parse_time(val):
    """Parse `val` into a duration in seconds.
    `val` can be a number or a string such as `'1y2m3d6h3m5s'`.
    Note that `'m'` is ambiguous between month and minute and will be interpreted depending on position (month by default)."""
    if isinstance(val, (float, int)):
        return val
    if not isinstance(val, str):
        raise ValueError('Time should be a number or a string.')
    regex = re.compile(r'((?P<y>\d+?)(Y|YR|YRS|YEARS|YEAR))?((?P<m>\d+?)(M|MONTHS|MONTH))?((?P<w>\d+?)(W|WEEKS|WEEK))?((?P<d>\d+?)D|DAY|DAYS)?((?P<h>\d+?)(H|HR|HRS|HOUR|HOURS))?((?P<min>\d+?)(M|MIN|MINUTE|MINUTES))?((?P<s>\d+?)(S|SEC|SECONDS|SECOND))?')
    val = val.replace(' ', '').upper()
    time = regex.match(val).groupdict()
    for key in time:
        if time[key] is None:
            time[key] = 0
        else:
            time[key] = float(time[key])
    secs = (((365 * time['y'] + 30 * time['m'] + 7 * time['w'] + time['d'])
            * 24 + time['h']) * 60 + time['min']) * 60 + time['s']
    return secs


def sec_to_str(val):
    """Return duration in seconds as a string."""
    duration = timedelta(seconds=val)
    return str(duration)


def generate_output_name(text='output'):
    """Generate output name based on current time."""
    time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    return f'{text}_{time}'


def check_path(path):
    """Check path and create missing folders if incorrect. Return path object."""
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)