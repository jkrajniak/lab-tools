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
import h5py

from md_libs import files_io


def _args():
    parser = argparse.ArgumentParser('Copy bond list from H5MD list to GROMACS topology')
    parser.add_argument('h5', help='Input H5MD file')
    parser.add_argument('topol', help='Input GROMACS Topology file')
    parser.add_argument('output', help='Output GROMACS topology file')
    parser.add_argument('--time_frame', help='Time frame', default=-1)
    parser.add_argument('--only_dynamic', action='store_true', default=False)

    return parser.parse_args()


def main():
    args = _args()

    h5 = h5py.File(args.h5, 'r')
    topol = files_io.GROMACSTopologyFile(args.topol)
    topol.read()

    if 'connectivity' in h5:
        for name, ds in list(h5['/connectivity'].items()):
            if not args.only_dynamic and isinstance(ds, h5py.Dataset):
                print(('Reading static {}'.format(name)))
                for b in ds:
                    if -1 not in b:
                        topol.new_data['bonds'][tuple(b)] = ['; h5md {}'.format(name)]
            elif isinstance(ds, h5py.Group):
                print(('Reading {}, time frame: {} of {}'.format(
                    name, 'last' if args.time_frame == -1 else args.time_frame, ds['step'].shape[0])))
                data = ds['value'][args.time_frame]
                for b in data:
                    if -1 not in b:
                        topol.new_data['bonds'][tuple(b)] = ['; h5md {}'.format(name)]

        topol.write(args.output)

if __name__ == '__main__':
    main()

