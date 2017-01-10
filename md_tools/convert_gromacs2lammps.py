#!/usr/bin/env python
"""
Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>

Backmapper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
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
import collections

from md_libs import files_io

__doc__ = 'Convert GROMACS to LAMMPS'


def prepare_bonded_lists(top):
    """Read GROMACS topology and group parameters into coefficient."""
    bonded_lists = {'bonds': top.bonds, 'angles': top.angles, 'dihedrals': top.dihedrals,
                    'improper_dihedrals': top.improper_dihedrals}
    return_bonded_lists = {'bonds': None, 'angles': None, 'dihedrals': None, 'improper_dihedrals': None}
    interaction_tuple = collections.namedtuple('InteractionBondedList', ['coeff', 'blist'])
    for bname, bond_list in bonded_lists.items():
        coeff = {}
        btypeid = 1
        particle_list = []
        for b, b_params in bond_list.items():
            if not b_params:
                raise RuntimeError('Empty bond params is not supported yet')
            key_params = tuple(
                b_params)  # tuple(['{:.3f}'.format(x) if x < 1.0 else '{:.1f}'.format(x) for x in map(float, b_params)])
            if coeff:
                interaction_id = coeff.setdefault(key_params, max(coeff.values())+1)
            else:
                interaction_id = coeff.setdefault(key_params, btypeid)
            particle_list.append(list(b) + [interaction_id])

        coeff = {tuple([int(float(k[0]))] + map(float, k[1:])): v for k, v in coeff.items()}

        return_bonded_lists[bname] = interaction_tuple(coeff, particle_list)

    return return_bonded_lists

def prepare_nonbonded_lists(top):
    """Read GROMACS topology and create atom types and mass lists."""

    def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

    atomtypename2typeid = {k['name']: i for i, k in enumerate(sorted(top.atomtypes.values()), 1)}
    top.atomtypename2typeid = atomtypename2typeid
    top.typeid2atomtypename = {}
    for type_name in top.atomtypes:
        top.atomtypes[type_name]['type_id'] = atomtypename2typeid[type_name]
        top.typeid2atomtypename[top.atomtypes[type_name]['type_id']] = type_name

    local_atomtypes = collections.defaultdict(dict)
    for at_id, at_obj in top.atoms.items():
        if at_obj.atom_type not in local_atomtypes:
            local_atomtypes[at_obj.atom_type] = {
                'mass': at_obj.mass,
                'charge': at_obj.charge,
                'sigma': 0.,
                'epsilon': 0.0}
        else:
            local_attype = local_atomtypes[at_obj.atom_type]
            if at_obj.mass != local_attype['mass']:
                print at_obj.mass, local_attype['mass']
    for at_type in local_atomtypes:
        if isclose(local_atomtypes[at_type]['mass'], top.atomtypes[at_type]['mass'], 0.0001, 0.0001):
            top.atomtypes[at_type]['mass'] = local_atomtypes[at_type]['mass']

    for at_id, at_obj in top.atoms.items():
        type_id = atomtypename2typeid[at_obj.atom_type]
        if not isclose(at_obj.mass, top.atomtypes[at_obj.atom_type]['mass'], 0.0001, 0.0001):
            new_type_id = max(atomtypename2typeid.values()) + 1
            top.atomtypes['{}{}'.format(at_obj.atom_type, new_type_id)] = top.atomtypes[at_obj.atom_type].copy()
            top.atomtypes['{}{}'.format(at_obj.atom_type, new_type_id)]['mass'] = at_obj.mass
            print('New atom type {}'.format(new_type_id))
            type_id = new_type_id
        at_obj.type_id = type_id


def prepare_atoms(in_top, in_gro):
    """Prepares atom list by reading coordinate and include in the topology object

    Args:
        in_top (files_io.GROMACSTopologyFile): input gromacs topology
    """
    for at_id in in_top.atoms:
        in_top.atoms[at_id].position = in_gro.atoms[at_id].position
    in_top.box = in_gro.box


def _args():
    parser = argparse.ArgumentParser(
        description='Convert GROMACS to LAMMPS topology',
        add_help=True)

    parser.add_argument('--in_gro', help='GRO file', required=True)
    parser.add_argument('--in_top', help='GROMACS top file', required=True)
    parser.add_argument('--out_lammps', help='Output LAMMPS data file')
    parser.add_argument('--scale_factor', default=10.0, help='Convert from nm to A', type=float)

    return parser


def main():
    args = _args().parse_args()

    in_gro = files_io.GROFile(args.in_gro)
    in_gro.read()
    in_top = files_io.GROMACSTopologyFile(args.in_top)
    in_top.read()
    in_top.replicate()

    bonded_lists = prepare_bonded_lists(in_top)
    for bname in bonded_lists:
        print('Found {} coeff of {} and {} interactions'.format(
            len(bonded_lists[bname].coeff.values()), bname, len(bonded_lists[bname].blist)))
    prepare_nonbonded_lists(in_top)
    prepare_atoms(in_top, in_gro)

    lmp_writer = files_io.LammpsWriter(args.out_lammps)
    lmp_writer.dist_scale = args.scale_factor
    lmp_writer.gro_topol = in_top
    lmp_writer.bonded_settings = bonded_lists
    lmp_writer.write()


if __name__ == '__main__':
    main()
