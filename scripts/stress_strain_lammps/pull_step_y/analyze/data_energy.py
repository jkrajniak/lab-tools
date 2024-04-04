#!/usr/bin/env python
"""
Copyright (C) 2017 Jakub Krajniak <jkrajniak@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import os
import sys
import multiprocessing as mp


def calculate_data(fname):
    columns = (0, 7, 8, 9, 10, 11)
    scale = np.array([1.0 for _ in columns])
    block_size = 100.0
    input_data = np.loadtxt(fname)[:, columns] * scale
    avg = np.average(input_data[:: int(block_size)], axis=0)
    std = np.std(input_data[:: int(block_size)], axis=0) * np.sqrt(1.0 / input_data[:: int(block_size)].shape[0])
    return [x for p in zip(avg, std) for x in p]


def sort_key(x):
    tmp = x.split(".")
    if len(tmp) == 4:
        if tmp[1] == "init":
            return -1
        else:
            return int(tmp[1])
    else:
        if tmp[2] == "init":
            return -1
        else:
            return int(tmp[2])


def main():
    xvg_files = [f for f in os.listdir(".") if f.endswith("xvg")]
    xvg_files.sort(key=sort_key)

    output = []

    p = mp.Pool()
    output = p.map(calculate_data, xvg_files)

    np.savetxt(sys.argv[1], output, header="s std pe pe_std ebond ebond_std eangle eangle_std edih edih_std eimp eimp_std")


if __name__ == "__main__":
    main()
