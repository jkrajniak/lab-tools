#!/usr/bin/env python
"""
Copyright (C) 2015 Jakub Krajniak <jkrajniak@gmail.com>

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
import h5py
import numpy as np

from matplotlib import pyplot as plt


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('h5file')
    parser.add_argument('--begin', default=0, type=int)
    parser.add_argument('--dt', default=1, type=float)

    return parser.parse_args()


def _print_observables(observables):
    for i, key in enumerate(observables, 1):
        print('{} - {}'.format(i, key))

    print('0 - Exit')


def _show_h5_observable(names, h5file, time_begin, dt, sum_values):
    plt.xlabel('Step')
    if sum_values:
        tvalues = []
        time = h5file['observables/{}/step'.format(names[0])][time_begin:]*dt
        for name in names:
            tvalues.append(h5file['observables/{}/value'.format(name)][time_begin:])
        values = np.sum(tvalues, axis=1)
        plt.plot(time, values)
    else:
        for name in names:
            values = h5file['observables/{}/value'.format(name)][time_begin:]
            time = h5file['observables/{}/step'.format(name)][time_begin:]*dt
            avg = np.average(values)

            print('Avg of {}: {}'.format(name, avg))

            plt.plot(time, values, label=name)

    plt.legend()
    plt.show()


def _show_csv_observable(names, data_file, time_begin, dt, sum_values, observables):

    plt.xlabel('Step')
    steps = data_file[:,0][time_begin:]*dt

    if sum_values:
        col_indexes = [observables.index(name) for name in names]
        values = np.sum(data_file[time_begin:, col_indexes], axis=1)
        plt.plot(steps, values)
    else:
        for name in names:
            col_index = observables.index(name)
            values = data_file[:,col_index][time_begin:]
            avg = np.average(values)

            print('Avg of {}: {}'.format(name, avg))

            plt.plot(steps, values, label=name)

    plt.legend()
    plt.show()


def main():
    args = _args()
    
    if args.h5file.endswith('h5'):
        data_file = h5py.File(args.h5file, 'r')
        observables = h5file['observables'].keys()
        _show_observable = _show_h5_observable
    else:
        data_file = np.loadtxt(args.h5file, skiprows=1)
        observables = open(args.h5file, 'r').readline().split()
        _show_observable = lambda w, x, y, z, s: _show_csv_observable(w, x, y, z, s, observables)
    
    _print_observables(observables)
    ans = raw_input('Select: ')
    while ans != '0':
        sum_values = False
        try:
            if '+' in ans:
                ans_index = map(int, ans.split('+'))
                sum_values = True
            else:
                ans_index = map(int, ans.split())
        except ValueError:
            ans = raw_input('Select: ')
            continue
        print sum_values
        obs_names = set([observables[x-1] for x in ans_index])
        _show_observable(obs_names, data_file, args.begin, args.dt, sum_values)
        ans = raw_input('Select: ')


if __name__ == '__main__':
    main()
