#!/usr/bin/env python
"""
Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>

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
import numpy as np

from matplotlib import pyplot as plt


def _args():
    parser = argparse.ArgumentParser('Calculate block average.')
    parser.add_argument('in_file', help='Input file, single column of time-series')
    parser.add_argument('--plot', '-p', help='Plot variance against block size',
                        action='store_true')
    parser.add_argument('--max_t', help='Maximum distance', type=int)
    parser.add_argument('--block', '-s', help='Block average size', type=float)

    return parser.parse_args()


def plot_blocks(data, max_t):
    mean = np.mean(data)
    var_tot = np.var(data)
    n = len(data)

    si = []
    for t in range(2, max_t):
        tb = np.floor(n/float(t))
        x_blocks = np.array_split(data, tb)
        v2 = np.average(np.power(list(map(np.mean, x_blocks)) - mean, 2))
        si.append([t, t*v2/var_tot])
    si = np.array(si)

    plt.plot(si[:, 0], si[:, 1])
    plt.xlabel('s')
    plt.ylabel(r'$s\sigma^2(<A_b>)/\sigma^2(A)$')
    plt.show()


def main():
    args = _args()

    data = np.loadtxt(args.in_file)

    if args.plot:
        plot_blocks(data, args.max_t)
    elif args.block:
        mean_value = np.mean(data)
        std_value = np.std(data, ddof=1)
        print(('AVG={}, std={}'.format(mean_value, std_value*np.sqrt(args.block/data.shape[0]))))


if __name__ == '__main__':
    main()
