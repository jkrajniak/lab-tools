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
import md_libs
import numpy as np
import collections


def _args():
    parser = argparse.ArgumentParser('List bonds from LAMMPS of given type')
    parser.add_argument('input_data', help='Input data file')
    parser.add_argument('out_list', help='Output file')
    parser.add_argument('--atom_type', help='Atom type',
                        required=True, type=int)
    parser.add_argument('--interactive', action='store_true')


    return parser.parse_args()


def main():
    args = _args()

    lammps_data = md_libs.files_io.LammpsReader()
    lammps_data.read_data(args.input_data)

    out_tuples = set()

    at_with_type = {k for k, v in lammps_data.atoms.items() if v['atom_type'] == args.atom_type}
    out_bonds = collections.defaultdict(list)
    for bl in lammps_data.topology['bonds'].values():
        for b1, b2 in bl:
            if b1 in at_with_type:
                out_bonds[b1].append(b2)
            if b2 in at_with_type:
                out_bonds[b2].append(b1)
            if b1 in at_with_type or b2 in at_with_type:
                out_tuples.add(tuple(sorted([b1, b2])))

    out_data = np.array(sorted(out_tuples))

    out_file = open(args.out_list, 'w')
    print('Write to {}'.format(args.out_list))
    out_file.writelines('\n'.join(['{} {}'.format(*x) for x in out_bonds.values()]))
    out_file.close()

    if args.interactive:
        import IPython;
        IPython.embed()


if __name__ == '__main__':
    main()

