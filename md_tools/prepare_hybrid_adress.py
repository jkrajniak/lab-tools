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
import collections
import copy
import numpy as np
import xml.etree.ElementTree as etree

from md_libs import files_io


BeadID = collections.namedtuple('BeadID', ['name', 'degree'])

XMLMap = collections.namedtuple(
    'XMLMap',
    ['ident',  # Identity name of molecules.
     'name',   # Name of molecules.
     'mass_map',  # Map with mass of CG beads.
     'molecule_beads',  # Returns molecule CG beads.
     'atom2cg',  # Map AT particle to CG bead
     'atom2cg_map',
     'at_topol',
     'cg_topol'
     ])

class BackmapperSettings(object):
    def __init__(self, input_xml):
        tree = etree.parse(input_xml)
        self.root = tree.getroot()
        self._parse()

    def _parse(self):
        cg_molecules = [self._parse_cg_molecule(r) for r in self.root.findall('cg_molecule')]
        self.cg_molecules = {r.name: r for r in cg_molecules}

    def _parse_cg_molecule(self, root):
        chain_name = root.find('name').text.strip()
        ident_name = root.find('ident').text.strip()

        at_topol = files_io.GROMACSTopologyFile(root.find('atomistic_topology').text.strip())
        at_topol.read()
        cg_topol = files_io.GROMACSTopologyFile(root.find('cg_topology').text.strip())
        cg_topol.read()

        # For each of cg bead and different degree holds a list of weights in the same
        # order as the <beads> section
        mapping2mass = {
            x.find('name').text: map(float, x.find('weights').text.strip().split())
            for x in root.iter('map')
        }
        # Prepares cg_bead_mass and cg_beads structures.
        cg_bead_mass = {}
        atom2cg_bead = {}
        cg_beads = {}
        atom2cg_map = {}
        for x in root.iter('cg_bead'):
            name = x.find('name').text.strip()
            atom_type = x.find('type').text.strip()
            try:
                degree = int(x.find('degree').text.strip())
            except:
                degree = '*'
            chain_name = chain_name
            mapping = x.find('mapping').text.strip()
            beads = x.find('beads').text.strip().split()
            active_site = x.find('active_site')
            if active_site is not None:
                active_site = active_site.text.strip()
            cg_bead_mass[BeadID(name, degree)] = {}
            atom_map_key = []
            for i, b_name in enumerate(beads):
                try:
                    cg_bead_mass[BeadID(name, degree)][b_name] = mapping2mass[mapping][i]
                    atom2cg_bead[b_name.split(':')[2]] = BeadID(name, degree)
                    atom_map_key.append(b_name.split(':')[2])
                except Exception as ex:
                    print i, b_name
                    raise ex
            atom2cg_map[tuple(atom_map_key)] = BeadID(name, degree)
            cg_beads[BeadID(name, degree)] = files_io.TopoAtom(
                atom_type=atom_type,
                name=name,
                chain_name=chain_name,
                active_site=active_site,
                mass=sum(cg_bead_mass[(name, degree)].values()))

        return XMLMap(
            ident=ident_name,
            name=chain_name,
            mass_map=cg_bead_mass,
            molecule_beads=cg_beads,
            atom2cg=atom2cg_bead,
            atom2cg_map=atom2cg_map,
            at_topol=at_topol,
            cg_topol=cg_topol)

def _args():
    parser = argparse.ArgumentParser('Creates AdResS coordinate file')
    parser.add_argument('--conf', required=True)
    parser.add_argument('--out', default='output.gro')
    parser.add_argument('--out_topol', default='output.top')
    parser.add_argument('--mapping', required=True)
    parser.add_argument('--interactive', action='store_true')

    return parser.parse_args()


def prepare_hybrid(settings, at_coordinate, args):
    """Prepare hybrid files."""
    out_coordinate = files_io.GROFile(args.out)
    out_coordinate.box = at_coordinate.box
    tmp_cg_mol = []
    tmp_cg_key = []
    new_at_id = 1
    cg_old2new_id = {}
    cg_new2old_id = {}
    at_old2new_id = {}
    at_new2old_id = {}
    at_cg_id = {}
    chain_name = None
    s = None
    for at in sorted(at_coordinate.atoms.keys()):
        at_data = at_coordinate.atoms[at]
        if chain_name is None:
            chain_name = at_data.chain_name
            s = settings.cg_molecules[chain_name]
        if chain_name != at_data.chain_name:
            raise RuntimeError('Only single chain name is supported, found {}'.format(at_data.chain_name))

        tmp_cg_key.append(at_data.name)
        tmp_cg_mol.append(at_data)

        if tuple(tmp_cg_key) in s.atom2cg_map:
            cg_bead = s.atom2cg_map[tuple(tmp_cg_key)]
            cg_at = s.cg_topol.chain_atom_names[chain_name][cg_bead.name][0]
            com = np.zeros(3)
            tot_mass = 0.0
            for a in tmp_cg_mol:
                com += s.at_topol.chain_atom_names[a.chain_name][a.name][0].mass * a.position
                tot_mass += s.at_topol.chain_atom_names[a.chain_name][a.name][0].mass
            com /= tot_mass
            out_coordinate.atoms[new_at_id] = files_io.Atom(
                atom_id=new_at_id,
                name=cg_at.name,
                chain_name=cg_at.chain_name,
                chain_idx=at_data.chain_idx,
                position=com
            )
            cg_id = new_at_id
            if at_data.chain_idx == 1:
                cg_old2new_id[cg_at.atom_id] = new_at_id
                cg_new2old_id[new_at_id] = cg_at.atom_id
            new_at_id += 1
            for a in tmp_cg_mol:
                out_coordinate.atoms[new_at_id] = files_io.Atom(
                            atom_id=new_at_id,
                            name=a.name,
                            chain_name=a.chain_name,
                            chain_idx=a.chain_idx,
                            position=a.position
                        )
                if a.chain_idx == 1:
                    at_old2new_id[a.atom_id] = new_at_id
                    at_new2old_id[new_at_id] = a.atom_id
                at_cg_id[new_at_id] = cg_id
                new_at_id += 1
            tmp_cg_key = []
            tmp_cg_mol = []

    out_coordinate.write(args.out, force=True)

    # Update topology
    output_topology = files_io.GROMACSTopologyFile(args.out_topol)
    output_topology.init(init_cross=True)

    for atid, at_data in out_coordinate.atoms.items():
        if at_data.chain_idx > 1:
            break
        s = settings.cg_molecules[at_data.chain_name]
        if atid in cg_new2old_id:
            cg_at = s.cg_topol.atoms[cg_new2old_id[atid]]
            output_topology.atoms[atid] = copy.copy(cg_at)
            output_topology.atoms[atid].atom_id = atid
        elif atid in at_new2old_id:
            at_at = s.at_topol.atoms[at_new2old_id[atid]]
            output_topology.atoms[atid] = copy.copy(at_at)
            output_topology.atoms[atid].atom_id = atid

    # Copy bonded sections
    s = settings.cg_molecules[chain_name]
    cg_topol = s.cg_topol
    output_topology.atomtypes = cg_topol.atomtypes
    for v in output_topology.atomtypes.values():
        v['type'] = 'V'

    output_topology.molecules = cg_topol.molecules
    output_topology.moleculetype = cg_topol.moleculetype

    for t, bparams in cg_topol.bonds.items():
        output_topology.cross_bonds[tuple(map(cg_old2new_id.get, t))] = bparams
    for t, aparams in cg_topol.angles.items():
        output_topology.cross_angles[tuple(map(cg_old2new_id.get, t))] = aparams
    for t, dparams in cg_topol.dihedrals.items():
        output_topology.cross_dihedrals[tuple(map(cg_old2new_id.get, t))] = dparams
    for t, pparams in cg_topol.pairs.items():
        output_topology.cross_dihedrals[tuple(map(cg_old2new_id.get, t))] = pparams

    # Update atomistic topology
    at_topol = s.at_topol
    def update_top_term(input_term, name):
        for t, bparams in input_term.items():
            new_t = tuple(map(at_old2new_id.get, t))
            cgs = map(at_cg_id.get, new_t)
            if cgs.count(cgs[0]) != len(cgs):  # Two different CG beads, cross bond
                output_topology.new_data['cross_{}'.format(name)][new_t] = bparams + [ '; cross at']
            else:
                output_topology.new_data[name][new_t] = bparams
    update_top_term(at_topol.bonds, 'bonds')
    update_top_term(at_topol.angles, 'angles')
    update_top_term(at_topol.dihedrals, 'dihedrals')
    update_top_term(at_topol.pairs, 'pairs')

    output_topology.write(args.out_topol)


def main():
    args = _args()

    at_coordinate = files_io.GROFile(args.conf)
    at_coordinate.read()

    settings = BackmapperSettings(args.mapping)
    prepare_hybrid(settings, at_coordinate, args)

    if args.interactive:
        import IPython
        IPython.embed()


if __name__ == '__main__':
    main()