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
import functools
import h5py
import numpy as np
import sys

from multiprocessing import Pool


from md_libs import _rdf


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('h5')
    parser.add_argument('--bins', default=100, type=int)
    parser.add_argument('-n', default=None, help='Index file')
    parser.add_argument('--type1', default=None, help='Particle type group 1', type=str)
    parser.add_argument('--type2', default=None, help='Particle type group 2', type=str)
    parser.add_argument('-b', default=0, type=int)
    parser.add_argument('-e', default=-1, type=int)
    parser.add_argument('--cutoff', default=-1, type=float)
    parser.add_argument('--no_normalize', default=False, action='store_true')
    parser.add_argument('--plot', action='store_true', default=False)
    parser.add_argument('--output', default=None, help='Output file')
    parser.add_argument('--nt', default=4, type=int)

    return parser.parse_args()


def get_single_rdf(h5file, type1, type2, index_file, cutoff, bins, do_norm, frame):
    h5 = h5py.File(h5file, 'r', driver='stdio', libver='latest')

    pids = None
    npart = None

    if index_file and (type1 or type2):
        print('Use index file or particle types, not both')
        sys.exit(1)

    has_types = False
    if index_file:
        with open(index_file, 'r') as findex:
            pids = list(map(int, ' '.join(findex.readlines()).split()))
            npart = len(pids)
    elif type1 is not None:
        has_types = True

    pos = h5['/particles/atoms/position/value']
    ids = h5['/particles/atoms/id/value']
    species = h5['/particles/atoms/species/value']
    L = h5['/particles/atoms/box/edges']
    if 'value' in L:
        L = L['value'][-1]
    vol = L[0] * L[1] * L[2]

    result = np.zeros(bins)
    multi = False
    print(frame)
    id_frame = ids[frame]
    p = pos[frame]
    npart = 0
    npart1 = 1
    npart2 = 1
    if has_types:
        species_frame = species[frame]
        pid_species = set()
        for t1 in type1:
            tt = id_frame[np.where(species_frame == t1)]
            pid_species.update(set(tt))
        pid_species = list(pid_species)
        p_pids = np.where(np.in1d(id_frame, pid_species))
        pp1 = p[p_pids]
        npart1 = len(set(pid_species))
        if type2 is not None:
            pid2_species = set()
            for t2 in type2:
                tt = id_frame[np.where(species_frame == t2)]
                pid2_species.update(set(tt))
            pid2_species = list(pid2_species)
            if set(pid2_species).intersection(set(pid_species)):
                raise RuntimeError('Particle ids of type1 and type2 overlap')
            p_pids2 = np.where(np.in1d(id_frame, pid2_species))
            pp2 = p[p_pids2]
            multi = True
            npart = len(set.union(set(pid_species), set(pid2_species)))
            npart2 = len(set(pid2_species))
        else:
            pp = pp1
            npart = len(set(pid_species))
    elif pids:
        p_pids = np.where(np.in1id(id_frame, pids))
        pp = p[p_pids]
    else:
        pp = p[np.where(id_frame != -1)]
        npart = len(set(id_frame[id_frame != -1]))
    if npart == 0:
        return None, None
    if multi:
        dx, tmp_r = _rdf.compute_rdf(
            np.asarray(pp1, dtype=np.double), np.asarray(pp2, dtype=np.double),
            L, bins, cutoff, False)
    else:
        dx, tmp_r = _rdf.compute_rdf(
            np.asarray(pp, dtype=np.double), None,
            L, bins, cutoff, False)
        npart2 = npart1

    phi = npart2/vol
    norm = phi*dx*npart1
    if do_norm:
        tmp_r /= norm
    result = np.nan_to_num(tmp_r)

    return result, dx*(np.arange(0, bins)+0.5)


def gets_rdf(h5, type1, type2, index_file, cutoff, bins=100, begin=0, end=-1, nt=4, do_norm=True):
    pos = h5['/particles/atoms/position/value']

    print(pos.shape[0])
    frames = list(range(begin, pos.shape[0] if end == -1 else end))

    if nt > 1:
        p = Pool(nt)

    get_rdf = functools.partial(get_single_rdf, h5.filename, type1, type2, index_file, cutoff, bins, do_norm)
    if nt > 1:
        results = p.map(get_rdf, frames)
    else:
        results = list(map(get_rdf, frames))

    x = results[0][1]
    result = np.array(np.sum([np.nan_to_num(k[0]) for k in results], axis=0))
    # This 
    norm = float(len(frames))

    return result/norm, x


def main():
    args = _args()
    h5 = h5py.File(args.h5, 'r')

    pids = None
    npart = None

    if args.n and (args.type1 or args.type2):
        print('Use index file or particle types, not both')
        sys.exit(1)

    has_types = False
    if args.n:
        with open(args.n, 'r') as findex:
            pids = list(map(int, ' '.join(findex.readlines()).split()))
            npart = len(pids)
    elif args.type1 is not None:
        has_types = True
        type1 = list(map(int, args.type1.split(',')))
        type2 = None
        if args.type2:
            type2 = list(map(int, args.type2.split(',')))

    result, x = gets_rdf(h5, type1, type2, args.n, args.cutoff, args.bins, args.b, args.e, args.nt, not args.no_normalize)
    if args.plot:
        from matplotlib import pyplot as plt
        plt.plot(x, result)
        plt.show()
    if args.output:
        np.savetxt(args.output, np.column_stack([x, result]))
        print(('File saved {}'.format(args.output)))


if __name__ == '__main__':
    main()
