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
        'This tool tries to make the particle index'
        ' continuous in the topology and in the coordinate files'))
    parser.add_argument('in_file', help='Input topology to fix')
    parser.add_argument('out_file', help='Input coordinate file')
    parser.add_argument('--remap_names', default=None)

    return parser.parse_args()


def replicate_list(input_list, id_map):
    return {tuple(map(id_map.get, k)): v for k, v in input_list.items()}


def fix_topol(in_topol, remap_names=None):
    if remap_names is None:
        remap_names = {}
    old2new_id = {}
    at_idx = 1
    new_atoms = {}
    last_chidx = -1
    ch_idx = 0
    for at_id in sorted(in_topol.atoms):
        data = in_topol.atoms[at_id]
        if last_chidx != data.chain_idx:
            ch_idx += 1
            last_chidx = data.chain_idx
        data.chain_idx = ch_idx
        old2new_id[data.atom_id] = at_idx
        data.atom_id = at_idx
        data.atom_name = remap_names.get(data.atom_name, data.atom_name)
        new_atoms[at_idx] = data
        at_idx += 1

    # Replicate lists.
    in_topol.bonds = replicate_list(in_topol.bonds, old2new_id)
    in_topol.angles = replicate_list(in_topol.angles, old2new_id)
    in_topol.dihedrals = replicate_list(in_topol.dihedrals, old2new_id)
    in_topol.improper_dihedrals = replicate_list(in_topol.improper_dihedrals, old2new_id)
    in_topol.pairs = replicate_list(in_topol.pairs, old2new_id)
    print([(k, v) for k, v in old2new_id.items() if k != v])


def fix_gro(in_gro, remap_names=None):
    if remap_names is None:
        remap_names = {}
    at_idx = 1
    new_atoms = {}
    last_chidx = -1
    ch_idx = 0
    for at_id in sorted(in_gro.atoms):
        if last_chidx != in_gro.atoms[at_id].chain_idx:
            ch_idx += 1
            last_chidx = in_gro.atoms[at_id].chain_idx
        at_data = in_gro.atoms[at_id]
        new_atoms[at_idx] = in_gro.atoms[at_id]._replace(atom_id=at_idx, chain_idx=ch_idx, name=remap_names.get(at_data.name, at_data.name))
        at_idx += 1
    in_gro.atoms = new_atoms

def main():
    args = _args()
    in_file = None
    topol = False
    remap_names = {}
    if args.remap_names:
        remap_names = dict(tuple(x.split(':')) for x in args.remap_names.split(','))
    if args.in_file.endswith('top'):
        in_file = files_io.GROMACSTopologyFile(args.in_file)
        in_file.read()
        topol = True
        fix_topol(in_file, remap_names)
        in_file.write(args.out_file, force=True)
    elif args.in_file.endswith('gro'):
        in_file = files_io.GROFile(args.in_file)
        in_file.read()
        fix_gro(in_file, remap_names)
        in_file.write(args.out_file, force=True)
    else:
        raise RuntimeError('Unknown input file')

if __name__ == '__main__':
    main()

