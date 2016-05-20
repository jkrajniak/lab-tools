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
import h5py
import itertools


def _args():
    parser = argparse.ArgumentParser('Shrink h5md file to given number of atoms')
    parser.add_argument('h5', help='Input h5md file')
    parser.add_argument('out_h5', help='Output file')
    parser.add_argument('--group', help='Atom group name', default='atoms')
    parser.add_argument('--N', help='Number of atoms (or -1 to auto value)', default=-1, type=int)

    return parser.parse_args()


def main():
    args = _args()
    in_h5 = h5py.File(args.h5, 'r')
    out_h5 = h5py.File(args.out_h5, 'w')

    for k in in_h5:
        in_h5.copy(k, out_h5)

    n_atoms = args.N
    if args.N == -1:
        # Determine number of atoms, mainly the 
        if 'id' not in in_h5['/particles/{}/'.format(args.group)]:
            raise RuntimeError('/id group not found, please define number of atoms')
        ids = in_h5['/particles/{}/id/value'.format(args.group)][-1]
        # Check if file is sorted
        ids_groups = [x[0] for x in itertools.groupby(ids)]
        if ids_groups.count(-1) > 1:
            raise RuntimeError('Is file sorted?')
        n_atoms = len([x for x in ids if x != -1])
    # Now clip datasets in /particles/<group>
    non_clip = ['box']
    for k, g in out_h5['/particles/{}'.format(args.group)].items():
        if k not in non_clip:
            print('Reshaping {}'.format(k))
            val = g['value']
            if len(val.shape) == 3:
                val.resize((val.shape[0], n_atoms, val.shape[2]))
            elif len(val.shape) == 2:
                val.resize((val.shape[0], n_atoms))
    out_h5.close()
    print('Saved in {}'.format(args.out_h5))

if __name__ == '__main__':
    main()

