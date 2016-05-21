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

parser = argparse.ArgumentParser(description='translate topology')
parser.add_argument('-t', '--top', help='GROMACS topology', required=True)
parser.add_argument('-i', '--index', help='Translate topology', default=0, type=int)
parser.add_argument('-o', '--output', help='Output', required=True)

args = parser.parse_args()

top = md_libs.files_io.GROMACSTopologyFile(args.top)
top.read()

max_atom_index = max(top.atoms)
min_atom_index = min(top.atoms)

replace_index = {k: k+args.index for k in range(min_atom_index, max_atom_index+1)}

output_file = open(args.output, 'w')

# Generate new atom entries
output_file.write('[ atoms ]\n')
for at_id in sorted(top.atoms):
  at = top.atoms[at_id]
  new_at_id = replace_index[at_id]
  fmt = '%s\t%s\t%s\t%s\t%s\t%s'
  data = [
      new_at_id,
      at.atom_type,
      at.chain_idx,
      at.chain_name,
      at.name,
      at.cgnr]
  if at.charge is not None:
    fmt += '\t%s'
    data.append(at.charge)
  if at.mass is not None:
    fmt += '\t%s'
    data.append(at.mass)
  fmt += '\n'
  output_file.write(fmt % tuple(data))

output_file.write('\n[ bonds ]\n')
new_bonds_dict = {
    tuple(map(replace_index.get, k)): v for k, v in top.bonds.iteritems()
    }
for k in sorted(new_bonds_dict):
  output_file.write('%s\t%s\n' % ('\t'.join(map(str, k)), ' '.join(new_bonds_dict[k])))

output_file.write('\n[ angles ]\n')
new_dict = {
    tuple(map(replace_index.get, k)): v for k, v in top.angles.iteritems()
    }
for k in sorted(new_dict):
  output_file.write('%s\t%s\n' % ('\t'.join(map(str, k)), ' '.join(new_dict[k])))

output_file.write('\n[ dihedrals ]\n')
new_dict = {
    tuple(map(replace_index.get, k)): v for k, v in top.dihedrals.iteritems()
    }
for k in sorted(new_dict):
  output_file.write('%s\t%s\n' % ('\t'.join(map(str, k)), ' '.join(new_dict[k])))

output_file.write('\n[ dihedrals ]\n')
new_dict = {
    tuple(map(replace_index.get, k)): v for k, v in top.improper_dihedrals.iteritems()
    }
for k in sorted(new_dict):
  output_file.write('%s\t%s\n' % ('\t'.join(map(str, k)), ' '.join(new_dict[k])))

output_file.write('\n[ pairs ]\n')
new_dict = {
    tuple(map(replace_index.get, k)): v for k, v in top.pairs.iteritems()
    }
for k in sorted(new_dict):
  output_file.write('%s\t%s\n' % ('\t'.join(map(str, k)), ' '.join(new_dict[k])))


output_file.close()
