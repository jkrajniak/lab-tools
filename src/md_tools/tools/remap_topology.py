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
import copy

from md_libs import files_io


def _args():
    parser = argparse.ArgumentParser('Remap topology partially with given id map')
    parser.add_argument('--in_top', required=True)
    parser.add_argument('--out_top', required=True)
    parser.add_argument('--mapping', required=True, help='The mapping file')

    return parser.parse_args()


def generate_new_bonded(in_bonded, mapping):
    output_data = {}
    for b_list, b_params in in_bonded:
        new_b_list = list(map(mapping.get, b_list))
        if new_b_list.count(None) == len(new_b_list):
            continue
        else:
            new_b_list = [mapping.get(p, p) for p in b_list]

        new_params = b_params[:]
        new_params.append('; get from {}'.format('-'.join(map(str, b_list))))
        output_data[tuple(new_b_list)] = new_params

    return output_data


def main():
    args = _args()

    in_top = files_io.GROMACSTopologyFile(args.in_top)
    in_top.read()

    mapping_file = open(args.mapping, 'r')

    mapping = {}
    new_names = {}
    for l in mapping_file.readlines():
        t = l.split()
        (old_id, new_id), new_name = list(map(int, t[:2])), t[2]
        if old_id not in in_top.atoms:
            raise RuntimeError('{} or {} not found in the atom list of topology'.format(
                old_id, new_id))
        if old_id in mapping:
            raise RuntimeError('Mapping file is wrong, {} already there'.format(old_id))
        mapping[old_id] = new_id
        new_names[new_id] = new_name
    print(('Read mapping file {}'.format(args.mapping)))

    # Prepare atoms if missing
    for old_id, new_id in sorted(list(mapping.items()), key=lambda x: x[1]):
        if new_id not in in_top.atoms:
            new_at = copy.copy(in_top.atoms[old_id])
            new_at.atom_id = new_id
            new_at.name = new_names[new_id]
            in_top.atoms[new_id] = new_at

    # Let's generate missing bonds, angles, dihedrals, pairs by taking
    in_top.bonds.update(generate_new_bonded(list(in_top.bonds.items()), mapping))
    in_top.angles.update(generate_new_bonded(list(in_top.angles.items()), mapping))
    in_top.dihedrals.update(generate_new_bonded(list(in_top.dihedrals.items()), mapping))
    in_top.improper_dihedrals.update(generate_new_bonded(list(in_top.improper_dihedrals.items()), mapping))
    in_top.pairs.update(generate_new_bonded(list(in_top.pairs.items()), mapping))

    in_top.write(args.out_top)

if __name__ == '__main__':
    main()

