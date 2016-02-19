#!/usr/bin/env python
"""
Copyright (C) 2015-2016 Jakub Krajniak <jkrajniak@gmail.com>

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
import cPickle
import h5py
import numpy as np
import sys

from libs import bonds
from libs import files_io


class ListAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        setattr(
            namespace,
            option_string.replace('-', ''),
            [map(int, x.split('-')) for x in values.split(',')]
        )


def _args():
    parser = argparse.ArgumentParser(('Calculating polymer distances, like'
                                      ' end-to-end distance or radius of gyration'))
    parser.add_argument('--trj', required=True, help='Input H5MD file')
    parser.add_argument('--molecules', help='Number of molecules', type=int, required=True)
    parser.add_argument('--N', type=int, help='Number of beads in molecule', required=True)
    parser.add_argument('--begin', '-b', type=int, default=0, help='Begin time frame')
    parser.add_argument('--end', '-e', type=int, default=-1, help='End time frame')
    parser.add_argument('--group', type=str, default='atoms',
                        help='Name of atom group')
    parser.add_argument('--out_ee', type=str, help='output end-to-end distance', default=None)
    parser.add_argument('--out_ee_points', type=str, help='Output end-to-end points', default=None)
    parser.add_argument('--out_rg', help='Output radius-of-gyration', default=None)
    parser.add_argument('--out_int', help='Output of average interal distance', default=None)
    return parser


def replicate_list(input_list, mol, N, shift=-1, cmplx=False):
    rr = [[i+x*N+shift for i in v] for x in range(mol) for v in input_list]
    if cmplx:
        return rr
    else:
        return [x for v in rr for x in v]


def calculate_rg(chains, chain_length, masses, box, half_box, tot_mass, frame):
    rg = np.zeros(chains)
    for ch in xrange(chains):  # Itrate through chains
        ch_com = np.sum(frame[ch*chain_length:ch*chain_length+chain_length], axis=0)
        ch_com /= float(chain_length)
        for n in xrange(chain_length):
            d = frame[ch*chain_length+n] - ch_com
            rg[ch] += masses[ch*chain_length+n]*np.sqrt(d.dot(d))
        rg[ch] /= tot_mass
    return rg


def fix_pbc(traj, box, half_box, chain_length, number_of_chains):
    print('Fix PBC for polymer chains')
    fidx = 0
    invBox = 1.0/box
    print invBox
    for frame in traj:
        sys.stdout.write('f={}\r'.format(fidx))
        sys.stdout.flush()
        for ch in xrange(number_of_chains):
            for n in xrange(chain_length-1):
                b1, b2 = frame[ch*chain_length+n], frame[ch*chain_length+n+1]
                for j in [0, 1, 2]:
                    d = b2[j] - b1[j]
                    b2[j] -= round(d*invBox[j])*box[j]
                d = b2 - b1
                if np.sqrt(d.dot(d)) > 0.5:
                    print ch*chain_length+n, ch*chain_length+n+1
                    sys.exit(1)
        fidx += 1
    return traj


def compute_ee(chains, chain_length, frame):
    ee = []
    for ch in xrange(chains):
        b1, b2 = frame[ch*chain_length], frame[ch*chain_length+chain_length-1]
        d = b2 - b1
        ee.append(np.sqrt(d.dot(d)))
    return ee


def main():
    args = _args().parse_args()

    data = h5py.File(args.trj, 'r')

    idx, box, trj, masses = files_io.prepare_h5md(data, args.group, args.begin, args.end,
                                                  no_image=True)

    half_box = 0.5*box
    trj = fix_pbc(trj, box, half_box, args.N, args.molecules)

    if args.out_ee:
        print('Calculating end-end distance...')
        ee = []
        for frame in trj:
            ee.extend(compute_ee(args.molecules, args.N, frame))
        out_ee = args.out_ee
        np.savetxt(out_ee, ee)
        print('Saved to {}'.format(out_ee))

    if args.out_ee_points:
        print('Save end-to-end points (raw data)...')
        ee_points = []
        for frame in trj:
            for ch in xrange(args.molecules):
                b1, b2 = frame[ch*args.N], frame[ch*args.N+args.N-1]
                ee_points.append([b1, b2])
        cPickle.dump(ee_points, open(args.out_ee_points, 'wb'))

    if args.out_rg:
        print('Calculating radius of gyration')
        rgs = []
        tot_masses = np.sum(masses[:args.N])
        for frame in trj:
            rg = calculate_rg(args.molecules, args.N, masses, box, half_box, tot_masses, frame)
            rgs.extend(rg)

        out_rg = args.out_rg
        np.savetxt(out_rg, rgs)
        print('Saved Rg^2 to {}'.format(out_rg))

    if args.out_int:
        print('Calculate internal distances.')
        int_distances = []
        fidx = 0
        for frame in trj:
            sys.stdout.write('{}\r'.format(fidx))
            sys.stdout.flush()
            int_distances.append(
                bonds.calculate_msd_internal_distance(frame, args.molecules, args.N, box, half_box))
            fidx += 1
        out_int = args.out_int
        np.savetxt(out_int, np.average(np.array(int_distances), axis=0))
        print('Saved internal distance to {}'.format(out_int))


if __name__ == '__main__':
    main()
