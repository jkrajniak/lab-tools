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

from md_libs import files_io


def _args():
    parser = argparse.ArgumentParser((
        'This tool tries to make the molecule index'
        ' continious in the topology and in the coordinate files'))
    parser.add_argument('in_file', help='Input topology to fix')
    parser.add_argument('out_file', help='Input coordinate file')

    return parser.parse_args()


def main():
    args = _args()
    in_file = None
    topol = False
    if args.in_file.endswith('top'):
        in_file = files_io.GROMACSTopologyFile(args.in_file)
        in_file.read()
        topol = True
    elif args.in_file.endswith('gro'):
        in_file = files_io.GROFile(args.in_file)
        in_file.read()
    else:
        raise RuntimeError('Unknown input file')

    last_chidx = -1
    ch_idx = 0
    for at_id in sorted(in_file.atoms):
        at_data = in_file.atoms[at_id]
        if last_chidx != at_data.chain_idx:
            ch_idx += 1
            last_chidx = at_data.chain_idx
        if topol:
            at_data.chain_idx = ch_idx
        else:
            in_file.atoms[at_id] = in_file.atoms[at_id]._replace(chain_idx=ch_idx)
        print in_file.atoms[at_id].chain_idx
    print('Saving updated file to {}'.format(args.out_file))
    in_file.write(args.out_file, force=True)

if __name__ == '__main__':
    main()

