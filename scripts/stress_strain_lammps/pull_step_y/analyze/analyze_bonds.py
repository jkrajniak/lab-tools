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
import numpy as np
import os
import cPickle
import collections
import multiprocessing as mp


def read_inputdata(input_file):
    timeseries_data = []
    with open(input_file, 'r') as inputf:
        timestep_item = False
        entries_item = False
        entries_header = None
        timestep_data = []
        for inline in inputf:
            if 'ITEM: TIMESTEP' in inline:
                timestep_item = True
                entries_item = False
                continue
            elif 'ITEM: ENTRIES' in inline:
                entries_item = True
                timestep_item = False
                entries_header = inline.split()[2:]
                continue
            if timestep_item:
                timestep = int(inline)
                timestep_item = False
                if not timestep_data:
                    continue
                timestep_data = np.array(timestep_data)
                tmp_data = np.zeros((timestep_data.shape[0], timestep_data.shape[1]+1))
                tmp_data[:, 0] = timestep
                tmp_data[:, 1:] = timestep_data
                timeseries_data.append(tmp_data)
                timestep_data = []
            elif entries_item:
                timestep_data.append(map(float, inline.split()[1:]))
    return timeseries_data


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('output_file')
    parser.add_argument('--single_file', default=None)
    parser.add_argument('--fileprefix', default='')
    parser.add_argument('--filesuffix', default='')
    parser.add_argument('--append', action='store_true')

    return parser.parse_args()


def process_file(filename):
    print('Processing {}'.format(filename))
    d = np.asfarray(read_inputdata(filename))
    if d.shape[0] == 0:
        print('File {} empty, skip'.format(filename))
        return (filename, None, None, None)
    bond_types = set(d[:,...,1].flat)
    print('Bond types: {}'.format(len(bond_types)))
    hist_per_bond = {}
    for bt in bond_types:
        hist_per_bond[bt] = np.histogram(d[d[:,...,1] == bt][:,...,4], bins='auto')
    global_hist = np.histogram(d[:,...,4], bins='auto')

    return (filename, bond_types, hist_per_bond, global_hist)



def main():
    args = _args()

    if args.single_file:
        output_data = [process_file(args.single_file)]
    else:
        pool = mp.Pool()
        output_data = pool.map(
            process_file,
            [f for f in os.listdir('.')
             if f.startswith(args.fileprefix) and f.endswith(args.filesuffix)])

    with open(args.output_file, 'wb') as output_file:
        cPickle.dump(output_data, output_file)

    print('Saved data to {}'.format(args.output_file))


if __name__ == '__main__':
    main()

