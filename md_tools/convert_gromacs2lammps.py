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
import datetime

from md_libs import files_io

__doc__ = 'Convert GROMACS to LAMMPS'


class LammpsWriter:
    def __init__(self, outfile):
        self.outfile = outfile
        self.gro_topol = None
        self.bonded_settings = None
        self.content = []
        self.dist_scale = 10.0  # nm -> A

    def kcal2kJ(self, v):
        """Converts kcal -> kJ"""
        return v * 4.184

    def kJ2kcal(self, v):
        """Converts kJ -> kcal"""
        return v / 4.184

    def _write_header(self):
        self.content.append(('LAMMPS data file via convert_gromacs2lammps, ver 1.0, {} '
                             ' topol file: {}\n').format(datetime.date.today(), self.gro_topol.file_name))
        self.content.append('{} atoms'.format(len(self.gro_topol.atoms)))
        self.content.append('{} atom types'.format(len(self.gro_topol.atomtypes)))
        self.content.append('{} bonds'.format(len(self.bonded_settings['bonds'].blist)))
        self.content.append('{} bond types'.format(len(self.bonded_settings['bonds'].coeff)))
        self.content.append('{} angles'.format(len(self.bonded_settings['angles'].blist)))
        self.content.append('{} angle types'.format(len(self.bonded_settings['angles'].coeff)))
        self.content.append('{} dihedrals'.format(len(self.bonded_settings['dihedrals'].blist)))
        self.content.append('{} dihedral types'.format(len(self.bonded_settings['dihedrals'].coeff)))
        self.content.append('{} impropers'.format(len(self.bonded_settings['improper_dihedrals'].blist)))
        self.content.append('{} improper types'.format(len(self.bonded_settings['improper_dihedrals'].coeff)))
        self.content.append('')

    def _write_box(self):
        box = self.gro_topol.box * self.dist_scale
        self.content.append('0 {} xlo xhi'.format(box[0]))
        self.content.append('0 {} ylo yhi'.format(box[1]))
        self.content.append('0 {} zlo zhi'.format(box[2]))
        self.content.append('')

    def _write_mass(self):
        self.content.append('Masses\n')
        for type_id in sorted(self.gro_topol.typeid2atomtypename):
            type_name = self.gro_topol.typeid2atomtypename[type_id]
            type_prop = self.gro_topol.atomtypes[type_name]
            self.content.append('{} {}'.format(type_id, type_prop['mass']))
        self.content.append('')

    def _write_pair_coeff(self):
        self.content.append('Pair Coeffs\n')
        for type_id in sorted(self.gro_topol.typeid2atomtypename):
            type_prop = self.gro_topol.atomtypes[self.gro_topol.typeid2atomtypename[type_id]]
            self.content.append('{} {} {}'.format(
                type_id, self.kJ2kcal(type_prop['epsilon']), type_prop['sigma']*self.dist_scale))
        self.content.append('')

    def _write_bond_coeffs(self):
        self.content.append('Bond Coeffs\n')
        btypeid_coeff = {v: k for k, v in self.bonded_settings['bonds'].coeff.items()}
        for btypeid in sorted(btypeid_coeff):
            coeff = btypeid_coeff[btypeid]
            if coeff[0] == 1:  # Harmonic
                K = 0.5*self.kJ2kcal(coeff[2])/(self.dist_scale**2)  # kJ/nm^2 -> kcal/A^2
                r0 = coeff[1]*self.dist_scale
                self.content.append('{} {:.3f} {:.3f}'.format(btypeid, K, r0))
            else:
                raise RuntimeError('Bond func type {} not supported yet'.format(coeff[0]))
        self.content.append('')

    def _write_angle_coeffs(self):
        self.content.append('Angle Coeffs\n')
        atypeid_coeff = {v: k for k, v in self.bonded_settings['angles'].coeff.items()}
        for atypeid in sorted(atypeid_coeff):
            coeff = atypeid_coeff[atypeid]
            if coeff[0] == 1:  # Harmonic
                K = 0.5*self.kJ2kcal(coeff[2])  # kJ/rad^2 -> kcal/rad^2
                theta0 = coeff[1]  # deg
                self.content.append('{} {:.3f} {:.3f}'.format(atypeid, K, theta0))
            else:
                raise RuntimeError('Angle func type {} not supported yet'.format(coeff[0]))
        self.content.append('')

    def _write_dihedral_coeffs(self):
        self.content.append('Dihedral Coeffs\n')
        dtypeid_coeff = {v: k for k, v in self.bonded_settings['dihedrals'].coeff.items()}
        for dtypeid in sorted(dtypeid_coeff):
            coeff = dtypeid_coeff[dtypeid]
            if coeff[0] == 3:  # RB -> nharmonic
                An = [pow(-1, n)*self.kJ2kcal(x) for n, x in enumerate(coeff[1:])]
                self.content.append('{} {} {}'.format(dtypeid, len(An), ' '.join(map('{:.3f}'.format, An))))
            else:
                raise RuntimeError('Dihedral func type {} not supported yet'.format(coeff[0]))
        self.content.append('')

    def _write_impropers_coeffs(self):
        self.content.append('Improper Coeffs\n')
        dtypeid_coeff = {v: k for k, v in self.bonded_settings['improper_dihedrals'].coeff.items()}
        for dtypeid in sorted(dtypeid_coeff):
            coeff = dtypeid_coeff[dtypeid]
            if coeff[0] == 1 and int(coeff[1]) == 180:  # harmonic improper to cvff
                K = self.kJ2kcal(coeff[2])
                d = -1
                n = int(coeff[3])
                self.content.append('{} {} {} {}'.format(dtypeid, K, d, n))
            else:
                raise RuntimeError('Improper func type {} not supported yet'.format(coeff[0]))
        self.content.append('')

    def _write_atoms(self):
        self.content.append('Atoms\n')
        last_res_id = 0
        res_id = 0
        for at_id in sorted(self.gro_topol.atoms):
            at_data = self.gro_topol.atoms[at_id]
            if at_data.chain_idx != last_res_id:
                res_id += 1
                last_res_id = at_data.chain_idx
            self.content.append('{at_id} {res_id} {type_id} {q} {x} {y} {z} 0 0 0'.format(
                at_id=at_id, res_id=res_id, type_id=at_data.type_id, q=at_data.charge,
                x=at_data.position[0]*self.dist_scale,
                y=at_data.position[1]*self.dist_scale,
                z=at_data.position[2]*self.dist_scale
            ))
        self.content.append('')

    def _write_bonded_list(self):
        bonded_list = [('Bonds', 'bonds'), ('Angles', 'angles'), ('Dihedrals', 'dihedrals'),
                       ('Impropers', 'improper_dihedrals')]
        for section_name, bonded_type in bonded_list:
            particle_list = sorted(self.bonded_settings[bonded_type].blist, key=lambda x: x[:-1])
            self.content.append('{}\n'.format(section_name))
            for idx, l in enumerate(particle_list, 1):
                type_id, plist = l[-1], l[:-1]
                self.content.append('{} {} {}'.format(idx, type_id, ' '.join(map(str, plist))))
            self.content.append('')


    def write(self):
        self._write_header()
        self._write_box()
        self._write_mass()
        self._write_pair_coeff()
        self._write_bond_coeffs()
        self._write_angle_coeffs()
        self._write_dihedral_coeffs()
        self._write_impropers_coeffs()
        self._write_atoms()
        self._write_bonded_list()

        with open(self.outfile, 'w') as outfile:
            outfile.write('\n'.join(self.content))


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
            key_params = tuple(['{:.3f}'.format(x) if x < 1.0 else '{:.1f}'.format(x) for x in map(float, b_params)])
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

    return parser


def main():
    args = _args().parse_args()

    in_gro = files_io.GROFile(args.in_gro)
    in_gro.read()
    in_top = files_io.GROMACSTopologyFile(args.in_top)
    in_top.read()

    bonded_lists = prepare_bonded_lists(in_top)
    for bname in bonded_lists:
        print('Found {} coeff of {} and {} interactions'.format(
            max(bonded_lists[bname].coeff.values()), bname, len(bonded_lists[bname].blist)))
    prepare_nonbonded_lists(in_top)
    prepare_atoms(in_top, in_gro)

    lmp_writer = LammpsWriter(args.out_lammps)
    lmp_writer.gro_topol = in_top
    lmp_writer.bonded_settings = bonded_lists
    lmp_writer.write()


if __name__ == '__main__':
    main()
