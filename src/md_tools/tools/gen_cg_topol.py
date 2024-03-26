#!/usr/bin/env python
"""
Copyright (C) 2017 Jakub Krajniak <jkrajniak@gmail.com>

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
import xml.etree.ElementTree as etree
from md_libs import files_io
import networkx as nx


class InputSettings:
    def __init__(self, input_file):
        self.input_file = input_file
        self.at2cg = {}
        self.mass_map = {}

    def parse(self):
        tree = etree.parse(self.input_file)
        root = tree.getroot()
        self.name = root.find('name').text
        for cg_bead in root.findall('*//cg_bead'):
            name = cg_bead.find('name').text
            beads = cg_bead.find('beads').text.split()
            cg_type = cg_bead.find('type').text
            for b in beads:
                _, chain_name, at_name = b.split(':')
                self.at2cg[at_name] = (name, cg_type)

        for mass_map in root.findall('*//map'):
            cg_type = mass_map.find('name').text
            weights = list(map(float, mass_map.find('weights').text.split()))
            self.mass_map[cg_type] = weights

def _args():
    parser = argparse.ArgumentParser('Map atomistic GROMACS topology to coarse-grained')
    parser.add_argument('--in_top', help='Input atomistic topology', required=True)
    parser.add_argument('--out_top', help='Output CG topology', required=True)
    parser.add_argument('--options', help='XML mapping files (comma separated) (VOTCA format)',
                        required=True)
    parser.add_argument('--itp', help='Optional ITP file with bondtypes, angletypes etc.',
                        required=False)

    return parser.parse_args()


def gen_angles_dihedrals(in_top, itp_top):
    g = nx.Graph()
    for b1, b2 in in_top.bonds:
        g.add_edge(b1, b2)

    dihedrals, angles = [], []
    for b1, b2 in g.edges():
        nbs1 = g.edge[b1]
        nbs2 = g.edge[b2]
        for nb2 in nbs2:
            if b1 == nb2:
                continue
            angles.append((b1, b2, nb2))
            nbs22 = g.edge[nb2]
            for nb22 in nbs22:
                if nb22 == b2:
                    continue
                dihedrals.append((b1, b2, nb2, nb22))
        for nb1 in nbs1:
            if nb1 == b2:
                continue
            angles.append((nb1, b1, b2))
            nbs11 = g.edge[nb1]
            for nb11 in nbs11:
                if nb11 == b1:
                    continue
                dihedrals.append((nb11, nb1, b1, b2))
        for nb1 in nbs1:
            for nb2 in nbs2:
                if nb1 == b2 or nb2 == b1:
                    continue
                dihedrals.append((nb1, b1, b2, nb2))
    # Unique
    angles = [a for a in angles[:] if list(reversed(a)) not in angles]
    angles.sort()
    # Add only internal angles.
    def add_values(input_list, param_list, output, only_internal):
        for k in input_list:
            rk = [in_top.atoms[p].chain_idx for p in k]
            if (rk.count(rk[0]) == len(rk)) == only_internal:
                rtypes = tuple([in_top.atoms[x].atom_type for x in k])
                try:
                    params = param_list[rtypes]
                    output[tuple(k)] = [params['func']] + params['params'] + ['; internal']
                except KeyError:
                    print(('Definition for {} not found'.format(rtypes)))
    add_values(angles, itp_top.angletypes, in_top.angles, True)
    add_values(angles, itp_top.angletypes, in_top.angles, False)

    dihedrals = [d for d in dihedrals[:] if list(reversed(d)) not in dihedrals]
    dihedrals.sort()
    add_values(dihedrals, itp_top.dihedraltypes, in_top.dihedrals, True)
    add_values(dihedrals, itp_top.dihedraltypes, in_top.dihedrals, False)

def map_topology(settings, in_top, itp_top):
    output_top = files_io.GROMACSTopologyFile('output')

    cg_beads, cg_bead = [], []
    current_cg_bead = None
    last_res_id, resid = 0, 1
    atid2res_id = {}
    molecules = {}
    for at_id in sorted(in_top.atoms):
        at_data = in_top.atoms[at_id]
        mol_settings = settings[at_data.chain_name]
        cg_name, cg_type = mol_settings.at2cg[at_data.name]
        cg_bead_mass = sum(mol_settings.mass_map[cg_type])
        if current_cg_bead is None:
            current_cg_bead = cg_name
            cg_bead = [cg_name, cg_type, cg_bead_mass,
                       resid, at_data.chain_idx, at_data.chain_name]
            last_res_id = at_data.chain_idx
        elif current_cg_bead != cg_name or last_res_id != at_data.chain_idx:
            cg_beads.append(cg_bead)
            resid += 1
            cg_bead = [cg_name, cg_type, cg_bead_mass,
                       resid, at_data.chain_idx, at_data.chain_name]
            current_cg_bead = cg_name
            last_res_id = at_data.chain_idx
        cg_bead.append(at_id)
        atid2res_id[at_id] = resid
    cg_beads.append(cg_bead)

    for cg_bead in cg_beads:
        cg_name, cg_type, cg_mass, cg_bead_id, res_id, res_name = cg_bead[:6]
        output_top.atoms[cg_bead_id] = files_io.TopoAtom(
            atom_id=cg_bead_id, atom_type=cg_type, chain_idx=res_id,
            chain_name=res_name, name=cg_name, cgnr=cg_bead_id,
            charge=0.0, mass=cg_mass)

    def gen_bond_list(input_list, output_list, input_types):
        for bs in sorted(input_list):
            res_ids = list(map(atid2res_id.get, bs))
            if res_ids.count(res_ids[0]) !=  len(res_ids):
                rtypes = tuple([output_top.atoms[x].atom_type for x in res_ids])
                chain_names = [output_top.atoms[x].chain_name for x in res_ids]
                try:
                    params = input_types[rtypes]
                except KeyError as ex:
                    print(('BS: {} res_ids: {} rtypes: {}'.format(
                        bs, res_ids, rtypes)))
                    raise ex
                bond_params = [params['func']] + params['params']
                bond_params.append(';')
                bond_params.append('-'.join(chain_names))
                if chain_names.count(chain_names[0]) != len(chain_names):
                    bond_params.append('cross-link')
                output_list[tuple(res_ids)] = bond_params

    gen_bond_list(in_top.bonds, output_top.bonds, itp_top.bondtypes)
    gen_angles_dihedrals(output_top, itp_top)

    return output_top

def main():
    args = _args()

    settings = {}
    for s in args.options.split(';'):
        st = InputSettings(s)
        st.parse()
        settings[st.name] = st

    in_top = files_io.GROMACSTopologyFile(args.in_top)
    in_top.read()

    itp_top = None
    if args.itp:
        itp_top = files_io.GROMACSTopologyFile(args.itp)
        itp_top.read()
    out_top = map_topology(settings, in_top, itp_top)
    out_top.write(args.out_top, force=True)

if __name__ == '__main__':
    main()

