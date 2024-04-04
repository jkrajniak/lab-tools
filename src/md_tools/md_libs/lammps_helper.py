"""
Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>

This file is part of lab-tools.

lab-tools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import StringIO
import collections
import math

import numpy as np


def get_lammps(filename, return_frames=False):
    """Gets data from LAMMPS log file

    Args:
        filename: The log filename.
        return_frames: Return the list of timeframes.

    Retunrs:
        The single StringIO object with the data or the list of StringIO objects if the data should be as the timeseries.

    """
    frame_lines = []
    with open(filename, "r") as of:
        lines = []
        in_section = False
        headers = []
        for l in of:
            if in_section and l.startswith("Loop"):
                in_section = False
                if return_frames:
                    frame_lines.append(lines)
                    lines = []
                continue
            if l.startswith("Time") or l.startswith("Step"):
                in_section = True
                headers.append(l.split())
                continue
            if in_section:
                lines.append(l)
        if return_frames:
            return [StringIO.StringIO("".join(l)) for l in frame_lines], headers
        else:
            return StringIO.StringIO("".join(lines)), headers[0]


def block_average(input_data, max_tb=200):
    """Calculated block average data.

    Args:
        input_data: The numpy array with the data to calculate.
        max_tb: Maximum block size

    Returns:
        The 2d numpy array, first column is the block size tb and second is the s.
    """
    n = input_data.shape[0]
    out = []

    tot_var = np.var(input_data)
    mean = np.mean(input_data)

    for tb in range(1, max_tb, 1):
        nb = int(math.ceil(n / tb))
        data = input_data  # [:tb*nb]
        nblocks = np.array_split(data, nb)
        var_block = np.mean([np.power(np.abs(np.average(b) - mean), 2) for b in nblocks])
        out.append([tb, tb * var_block / tot_var])
    out = np.array(out)
    return out


def parse_timedata(filename):
    """Parse time data from dump command."""
    timeframes = collections.defaultdict(list)
    with open(filename, "r") as fo:
        first_frame = True
        time_step = None
        nrows = 0
        for line in fo:
            if line.startswith("#"):
                continue
            if first_frame:
                time_step, nrows = map(int, line.split())
                first_frame = False
                rowcounter = 0
                continue
            if rowcounter < nrows:
                timeframes[int(time_step)].append(map(float, line.split()))
                rowcounter += 1
                continue
            if rowcounter == nrows:
                time_step, nrows = map(int, line.split())
                rowcounter = 0
    for time_frame in timeframes:
        timeframes[time_frame] = np.array(timeframes[time_frame])
    return timeframes


def load_numpyarray(filename):
    """Loads data from file intro struct numpy array."""
    data = np.loadtxt(filename, skiprows=1)
    header = [(x, np.float) for x in open(filename).readline().replace("# ", "").split()]
    data.dtype = header
    return data
