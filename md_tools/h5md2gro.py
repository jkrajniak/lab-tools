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

from md_libs import files_io


def _args():
    parser = argparse.ArgumentParser('Update atom positions of input gro file from H5MD file.')
    parser.add_argument('--h5', required=True)
    parser.add_argument('--group', required=True)
    parser.add_argument('--input_gro', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--frame', type=int, default=-1)
    parser.add_argument('--unfolded', action='store_true', default=False)
    parser.add_argument('--valid_species')
    parser.add_argument('--scale_factor', default=1.0, type=float)
    return parser.parse_args()


def main():
    args = _args()

    h5 = h5py.File(args.h5)
    in_gro = files_io.GROFile(args.input_gro)
    in_gro.scale_factor = args.scale_factor
    in_gro.read()

    pos = h5['/particles/{}/position/value'.format(args.group)][args.frame]
    try:
        species = h5['/particles/{}/species/value'.format(args.group)][args.frame]
    except:
        species = h5['/particles/{}/species'.format(args.group)][args.frame]

    valid_species = None
    if args.valid_species:
        valid_species = set(map(int, args.valid_species.split(',')))

    images = h5['/particles/{}/image/value'.format(args.group)][args.frame]
    print images
    try:
        box = np.array(h5['/particles/{}/box/edges/value'.format(args.group)][args.frame])
    except:
        box = np.array(h5['/particles/{}/box/edges'.format(args.group)][args.frame])
    print('Box: {}'.format(box))

    ids = sorted(in_gro.atoms)
    ppid = 0
    for pid, p in enumerate(pos):
        if pid >= len(in_gro.atoms):
            break
        if valid_species is None or species[pid] in valid_species:
            at_data = in_gro.atoms[ids[ppid]]
            if args.unfolded:
                in_gro.atoms[ids[ppid]] = at_data._replace(position = p + images[pid]*box)
            else:
                in_gro.atoms[ids[ppid]].position = at_data._replace(position=p)
            ppid += 1
    in_gro.write(args.output, force=True)


if __name__ == '__main__':
    main()
