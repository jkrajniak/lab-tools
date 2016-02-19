#!/usr/bin/env python
"""
Copyright (C) 2014 Jakub Krajniak <jkrajniak@gmail.com>

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



def _args():
    parser = argparse.ArgumentParser('Compute histogram and save to the disk')
    parser.add_argument('--bins', help='Number of bins', default=100, type=int)
    parser.add_argument('--normed', help='Normed', default=1, type=int)
    parser.add_argument('--plot', help='Plot?', action='store_true', default=False)
    parser.add_argument('--scale', help='Multiply values by', default=1.0, type=float)
    parser.add_argument('--interactive', help='Interactive', action='store_true', default=False)
    parser.add_argument('--type', choices=('bonds', 'angles', 'no'), default='no')
    parser.add_argument('input_file', help='Input file', nargs='+')
    parser.add_argument('output_file', help='Output file', default=None, nargs='?')
    return parser.parse_args()


def main():
    args = _args()
    raw_data = []
    for input_file in args.input_file:
        print('Reading file {}'.format(input_file))
        data = np.loadtxt(input_file) * args.scale
        raw_data.extend(data)
    raw_data = np.array(raw_data)

    data = raw_data[np.isfinite(raw_data)]

    if len(raw_data) != len(data):
        print('Warning! Raw data contains NaN and Infinite values ({})'.format(
            len(raw_data)-len(data)))

    n, bins = np.histogram(data, bins=args.bins, density=args.normed)

    bw = bins[2] - bins[1]
    bw = bins + bw
    if args.type == 'bonds':
        n /= (4*np.pi*(bw[:len(n)])**2)
    elif args.type == 'angles':
        n /= np.rad2deg(np.sin(np.deg2rad(bw[:len(n)])))
    dt = bins[2] - bins[1]
    n_sum = sum(n*(bins[2]-bins[1]))

    print('A bit of statistics\n====================')
    print('Number of samples {}'.format(len(data)))
    print('Avg: {}'.format(np.average(data)))
    print('Var: {}'.format(np.var(data)))
    print('Std: {}'.format(np.std(data)))
    print('Sum: {}'.format(n_sum))
    n_list = list(n)
    n_max = max(n_list)
    n_index = n_list.index(n_max)
    print('Max: {}'.format(bins[n_index]+dt))

    if args.interactive:
        import IPython
        IPython.embed()

    if args.plot:
	from matplotlib import pyplot as plt
        plt.plot([x+dt for x in bins][:-1], n)
        plt.show()

    if args.output_file is None:
        args.output_file = '{}.hist'.format(args.input_file[0].split('.')[:-1][0])

    print('Save output to {}'.format(args.output_file))
    np.savetxt(args.output_file, zip([x+dt for x in bins][:-1], n))


if __name__ == '__main__':
    main()
