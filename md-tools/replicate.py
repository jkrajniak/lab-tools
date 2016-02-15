"""
Copyright (C) 2014 Jakub Krajniak <jkrajniak@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


The idea of this program is to create topology file that will be ready to do the cross-link
operation. To achive this we need to explicitly define the bonds between the atoms.
"""

import argparse

from libs import files_io

print 'replicate'
print 'usage: '
print '  replicate.py -c MEL:Mela.itp,FOR:For.itp -f conf.gro -o output.top'

parser = argparse.ArgumentParser(description='Create topology for mixture')
parser.add_argument('-c', '--config',
    help='Define the config with the format chain_id:topology_file,chain_id:topology_file',
    required=True)
parser.add_argument('-f', '--pdb',
    help='The coordinate file with the mixture.',
    required=True
    )
parser.add_argument('-o', '--out',
    help='The output topology file.',
    required=True
    )


args = parser.parse_args()
config = args.config

molecules = config.split(',')
mol_top = {}
mol_atom_nr = {}
print 'Config:'
for mol in molecules:
  print '\t-', mol
  mol_name, top_file = mol.split(':')
  mol_top[mol_name] = files_io.TopologyFile(top_file)
  mol_top[mol_name].open()
  mol_top[mol_name].read()
  mol_atom_nr[mol_name] = len(mol_top[mol_name].atoms)

if args.pdb.endswith('pdb'):
  coordinate_file = files_io.PDBFile(args.pdb)
elif args.pdb.endswith('gro'):
  coordinate_file = files_io.GROFile(args.pdb)
else:
  raise Exception('Unsupported coordinate file format.')

coordinate_file.open()
coordinate_file.read()

# Output topology
output_top = files_io.TopologyFile(args.out)

# Fill up the new topology by reading the coordinate file and use the atom_id from that file.
# We assume that each of the topology forms consistent block of entries so there is no mix.
current_ch_idx, previous_ch_idx = None, None
for atom_id in sorted(coordinate_file.data):
  at = coordinate_file.data[atom_id]
  previous_ch_idx = current_ch_idx
  current_ch_idx = at.chain_idx
  if previous_ch_idx != current_ch_idx: # New chain
    top = mol_top[at.chain_name]
    # Replicate atoms
    new_at_id = at.atom_id
    at_id_map = {}  # The atom id map
    for at_id in top.atoms:
      atom_line = top.atoms[at_id].get_tuple()
      atom_line[0] = new_at_id
      atom_line[2] = current_ch_idx
      atom_line[5] = new_at_id
      output_top.new_data['atoms'][new_at_id] = atom_line
      at_id_map[at_id] = new_at_id
      new_at_id += 1

    # Replicate bonds
    for b, b_def in top.bonds.iteritems():
      b_tuple = tuple(map(at_id_map.get, b))
      output_top.new_data['bonds'][b_tuple] = b_def

    for a, a_def in top.angles.iteritems():
      a_tuple = tuple(map(at_id_map.get, a))
      output_top.new_data['angles'][a_tuple] = a_def

    for d, d_def in top.dihedrals.iteritems():
      d_tuple = tuple(map(at_id_map.get, d))
      output_top.new_data['dihedrals'][d_tuple] = d_def

    for p, p_def in top.pairs.iteritems():
      p_tuple = tuple(map(at_id_map.get, p))
      output_top.new_data['pairs'][p_tuple] = p_def

    for imp_d, imp_d_def in top.improper_dihedrals.iteritems():
      imp_d_tuple = tuple(map(at_id_map.get, imp_d))
      output_top.new_data['improper_dihedrals'][imp_d_tuple] = imp_d_def

output_top.write(args.out)
