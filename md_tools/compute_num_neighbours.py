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
import h5py
import itertools
import numpy as np
from scipy.integrate import quad

from md_libs import _rdf

from matplotlib import pyplot as plt


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('h5', help='Input HDF5 file')
    parser.add_argument('--cutoff', default=None, help='Cutoff distance', type=float)
    parser.add_argument('--begin', '-b', default=0, help='Begin frame', type=int)
    parser.add_argument('--end', '-e', default=-1, help='End frame', type=int)
    parser.add_argument('--type1', '-t1', type=str, required=True,
                        help='Types 1 (separated by comma)')
    parser.add_argument('--type2', '-t2', type=str, required=True,
                        help='Types 2 (separated by comma)')

    return parser.parse_args()


def main():
    args = _args()

    args = _args()
    h5 = h5py.File(args.h5, 'r')

    pos = h5['/particles/atoms/position/value']
    ids = h5['/particles/atoms/id/value']
    species = h5['/particles/atoms/species/value']
    #states = h5['/particles/atoms/state/value']

    L = h5['/particles/atoms/box/edges']
    if 'value' in L:
        L = L['value'][-1]

    vol = L[0]*L[1]*L[2]

    result = []
    for frame in xrange(args.begin, pos.shape[0] if args.end == -1 else args.end):
        id_frame = ids[frame]
        p = pos[frame]
        species_frame = species[frame]
        pid_species1 = set()
        for t1 in map(int, args.type1.split(',')):
            tt = id_frame[np.where(species_frame == t1)]
            pid_species1.update(set(tt))
        pid_species2 = set()
        for t2 in map(int, args.type2.split(',')):
            tt = id_frame[np.where(species_frame == t2)]
            pid_species2.update(set(tt))

        pp1 = p[np.in1d(id_frame, list(pid_species1))]
        pp2 = p[np.in1d(id_frame, list(pid_species2))]
        #pp1 = p[np.where(id_frame != -1)]
        #pp2 = p[np.where(id_frame != -1)]

        avg_num = _rdf.compute_nb(
            np.asarray(pp1, dtype=np.float),
            np.asarray(pp2, dtype=np.float), L, args.cutoff)

        result.append([frame, np.average(avg_num)])

    result = np.array(result)
    plt.plot(result[:, 0], result[:, 1])
    plt.show()
    print(np.average(result[:, 1]))

if __name__ == '__main__':
    main()

