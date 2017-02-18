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

from multiprocessing import Pool
import functools

#from matplotlib import pyplot as plt


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
    parser.add_argument('--output', required=True)
    parser.add_argument('--nt', type=int, default=4, help='Number of CPUs')

    return parser.parse_args()


def get_avg_nb(types1, types2, L, cutoff, ids, pos, species, frame):
    id_frame = ids[frame]
    p = pos[frame]
    species_frame = species[frame]
    pid_species1 = set()
    for t1 in types1:
        tt = id_frame[np.where(species_frame == t1)]
        pid_species1.update(set(tt))
    pid_species2 = set()
    for t2 in types2:
        tt = id_frame[np.where(species_frame == t2)]
        pid_species2.update(set(tt))

    pp1 = p[np.in1d(id_frame, list(pid_species1))]
    pp2 = p[np.in1d(id_frame, list(pid_species2))]

    avg_num = _rdf.compute_nb(
        np.asarray(pp1, dtype=np.float),
        np.asarray(pp2, dtype=np.float), L, cutoff)

    print frame
    return (frame, np.average(avg_num))



def get_avg_nb2(types1, types2, L, cutoff, filename, frame):
    h5 = h5py.File(filename, 'r', libver='latest', driver='stdio')

    pos = h5['/particles/atoms/position/value']
    ids = h5['/particles/atoms/id/value']
    species = h5['/particles/atoms/species/value']

    id_frame = ids[frame]
    p = pos[frame]
    species_frame = species[frame]
    pid_species1 = set()
    for t1 in types1:
        tt = id_frame[np.where(species_frame == t1)]
        pid_species1.update(set(tt))
    pid_species2 = set()
    for t2 in types2:
        tt = id_frame[np.where(species_frame == t2)]
        pid_species2.update(set(tt))

    pp1 = p[np.in1d(id_frame, list(pid_species1))]
    pp2 = p[np.in1d(id_frame, list(pid_species2))]

    avg_num = _rdf.compute_nb(
        np.asarray(pp1, dtype=np.float),
        np.asarray(pp2, dtype=np.float), L, cutoff)

    #h5.close()

    print frame
    return (frame, np.average(avg_num))


def main():
    args = _args()
    h5filename = args.h5

    h5 = h5py.File(h5filename, 'r')

    pos = h5['/particles/atoms/position/value']
    #ids = h5['/particles/atoms/id/value']
    #species = h5['/particles/atoms/species/value']
    #states = h5['/particles/atoms/state/value']

    L = h5['/particles/atoms/box/edges']
    if 'value' in L:
        L = L['value'][-1]

    L = np.asarray(L, dtype=np.double)

    cutoff = args.cutoff
    if args.cutoff is None:
        cutoff = 0.5*L[0]

    print('L: {}, cutoff: {}'.format(L, cutoff))

    types1 = map(int, args.type1.split(','))
    types2 = map(int, args.type2.split(','))

    p = Pool(args.nt)

    result = []

    frames = range(args.begin, pos.shape[0] if args.end == -1 else args.end)
    frames_split = np.array_split(frames, 100)
    h5.close()

    #for frames in frames_split:
    #    for frame in frames:
    #        result.append(get_avg_nb(types1, types2, L, cutoff, ids, pos, species, frame))

    #get_avg_nb_ = functools.partial(get_avg_nb, types1, types2, L, cutoff, np.asarray(ids), np.asarray(pos), np.asarray(species))
    get_avg_nb_ = functools.partial(get_avg_nb2, types1, types2, L, cutoff, h5filename)
    #print 'bcc'
    result = p.map(get_avg_nb_, frames)

    #for frame in xrange(args.begin, pos.shape[0] if args.end == -1 else args.end):
    #    result.append([frame, get_avg_nb(types1, types2, L, cutoff, ids, pos, species, frame)])
    #    print frame

    result = np.array(result)
    np.savetxt(args.output, result)
    print('Saved data {}'.format(args.output))
    print('Average neighbours {}'.format(np.average(result, axis=0)))

if __name__ == '__main__':
    main()
