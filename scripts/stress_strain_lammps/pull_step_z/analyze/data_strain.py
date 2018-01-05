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

import argparse
import math
import numpy as np
import os
import sys
import multiprocessing as mp


def calculate_data(fname):
    columns = (0, 1, 2, 3, 4, 5, 6)
    scale = np.array([1.0, -0.101325, -0.101325, -0.101325, 1.0, 1.0, 1.0])
    input_data = np.loadtxt(fname)[:, columns]*scale
    data_length = input_data.shape[0]
    input_data = input_data[int(data_length/2):]

    number_of_samples = 100
    block_size = math.ceil(input_data.shape[0]/number_of_samples)
    if block_size < 1.0:
        block_size = 1
    print('bloc size: {}'.format(block_size))
    print('data size: {}'.format(input_data.shape[0]))
    avg = np.average(input_data[::int(block_size)], axis=0)
    std = np.std(input_data[::int(block_size)], axis=0)/np.sqrt(input_data[::int(block_size)].shape[0]-1)
    return [x for p in zip(avg, std) for x in p]

def sort_key(x):
    tmp = x.split('.')
    if len(tmp) == 4:
        if tmp[1] == 'init':
            return -1
        else:
            return int(tmp[1])
    else:
        if tmp[2] == 'init':
            return -1
        else:
            return int(tmp[2])

def main():
    xvg_files = [f for f in os.listdir('.') if f.endswith('xvg')]

    xvg_files.sort(key=sort_key)

    output = []

    #p = mp.Pool()
    for xvg in xvg_files:
        print('File: {}'.format(xvg))
        output.append(calculate_data(xvg))
    #output = p.map(calculate_data, xvg_files)

    np.savetxt(sys.argv[1], output, header='s std pxx pxx_std pyy pyy_std pzz pzz_std lx lx_std ly ly_std lz lz_std')


if __name__ == '__main__':
    main()
