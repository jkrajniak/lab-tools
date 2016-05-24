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
from md_libs import files_io

__doc__ = 'Reorder atoms in topology, to match the order in coordinate file.'


def _args():
    parser = argparse.ArgumentParser('Reorder atoms in topology to match what is in coord file.')
    parser.add_argument('--in_top', required=True)
    parser.add_argument('--out_top', required=True)
    parser.add_argument('--coord', required=True)
    parser.add_argument('--clean', action='store_true', default=False)
    parser.add_argument('--remove_cross', action='store_true', default=False)

    return parser.parse_args()


def generate_list(input_list, id_map, clean=False):
    output = {}
    if clean:
        for k, v in input_list.items():
            try:
                new_k = tuple(map(lambda x: id_map[x], k))
                output[new_k] = v
            except KeyError:
                continue
    else:
        output = {tuple(map(lambda x: id_map[x], k)): v
                 for k, v in input_list.items()}
    return output

def main():
    args = _args()

    in_top = files_io.GROMACSTopologyFile(args.in_top)
    in_top.read()

    coord = files_io.GROFile(args.coord)
    coord.read()

    # chain_idx:chain_name:atom_name
    input_name2id = {
        '{}:{}:{}'.format(v.chain_idx, v.chain_name, v.name): k
        for k, v in coord.atoms.iteritems()}
    topol_name2id = {
        '{}:{}:{}'.format(v.chain_idx, v.chain_name, v.name): k
        for k, v in in_top.atoms.iteritems()}

    # Map topol id -> atom_id
    topol_old2new = {}
    if args.clean:
        for x in input_name2id:
            topol_old2new[topol_name2id[x]] = input_name2id[x]
    else:
        topol_old2new = {
            topol_name2id[x]: input_name2id[x] for x in topol_name2id
        }

    new_topol_atoms = {}
    for x in in_top.atoms:
        at = in_top.atoms[x]
        try:
            at.atom_id = topol_old2new[x]
        except KeyError:
            print('Skiping atom {}:{}'.format(at.chain_name, at.name))
            continue
        new_topol_atoms[topol_old2new[x]] = at

    in_top.atoms = new_topol_atoms

    new_bonds = generate_list(in_top.bonds, topol_old2new, args.clean)
    new_angles = generate_list(in_top.angles, topol_old2new, args.clean)
    new_dihs = generate_list(in_top.dihedrals, topol_old2new, args.clean)
    new_pairs = generate_list(in_top.pairs, topol_old2new, args.clean)
    new_cr_bonds = generate_list(in_top.cross_bonds, topol_old2new, args.clean)
    new_cr_angles = generate_list(in_top.cross_angles, topol_old2new, args.clean)
    new_cr_dihs = generate_list(in_top.cross_dihedrals, topol_old2new, args.clean)
    new_cr_pairs = generate_list(in_top.cross_pairs, topol_old2new, args.clean)

    in_top.bonds = new_bonds
    in_top.angles = new_angles
    in_top.dihedrals = new_dihs
    in_top.pairs = new_pairs

    if args.remove_cross:
        in_top.bonds.update(new_cr_bonds)
        in_top.angles.update(new_cr_angles)
        in_top.dihedrals.update(new_cr_dihs)
        in_top.pairs.update(new_cr_pairs)
        in_top.cross_bonds = {}
        in_top.cross_angles = {}
        in_top.cross_dihedrals = {}
        in_top.cross_pairs = {}
    else:
        in_top.cross_bonds = new_cr_bonds
        in_top.cross_angles = new_cr_angles
        in_top.cross_pairs = new_cr_pairs
        in_top.cross_dihedrals = new_cr_dihs

    in_top.write(args.out_top, force=True)

if __name__ == '__main__':
    main()
