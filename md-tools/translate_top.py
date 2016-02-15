
import argparse

from libs import files_io

parser = argparse.ArgumentParser(description='translate topology')
parser.add_argument('-t', '--top', help='GROMACS topology', required=True)
parser.add_argument('-i', '--index', help='Translate topology', default=0, type=int)
parser.add_argument('-o', '--output', help='Output', required=True)

args = parser.parse_args()

top = files_io.TopologyFile(args.top)
top.open()
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
