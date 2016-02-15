"""
Copyright (C) 2014 Jakub Krajniak <jkrajniak@gmail.com>

This file is part of BondMatcher.

BondMatcher is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

BondMatcher - I/O library. Handles opening and writing different files."""


from collections import defaultdict
import numpy
import os
import sys

import structures  # pylint:disable=W0403


def prepare_path(file_path):
  """Prepare the file to open.

  Args:
    file_path: The file path.

  Returns:
    The path to the file.
  """

  if os.path.exists(file_path):
    file_name = os.path.basename(file_path)
    dir_name = os.path.dirname(file_path)
    if not dir_name:
      dir_name = '.'
    existing_copies = [x for x in os.listdir(dir_name) if x.startswith('_%s' % file_name)]
    if existing_copies:
      max_copy_id = max([int(x.strip('_').split('.')[-1]) for x in existing_copies])
    else:
      max_copy_id = 0
    new_file_name = '_%s.%d_' % (file_name, max_copy_id+1)
    new_file_path = os.path.join(dir_name, new_file_name)
    print '\nFound: %s, backup on: %s\n' % (file_path, new_file_path)
    os.rename(file_path, new_file_path)

  return file_path


def sort_h5md_array(input_array, ids, max_T = None):
    """Sorts H5MD dataset"""

    T = len(input_array)
    output_shape = list(input_array.shape)
    if max_T:
        T = max_T
        output_shape[0] = T
    output_array = numpy.zeros(output_shape)
    # Iterates over time frames.
    for t in xrange(T):
        sys.stdout.write('Progress: {:.2f} %\r'.format(100.0*float(t)/T))
        sys.stdout.flush()
        idd = [
	    x[1] for x in sorted(
	        [(p_id, col_id) for col_id, p_id in enumerate(ids[t])],
	        key=lambda y: (True, y[0]) if y[0] == -1 else (False, y[0]))
	    ]
        output_array[t] = input_array[t][idd]
    return output_array

def prepare_h5md(h5file, group_name, begin, end, step=None, no_image=False):
    """Returns H5MD data that are sorted and transformed."""

    if step is None:
        step = 1

    # Checks if there is an ids data set. That implies sorting.
    ids = None
    if 'id' in h5file['/particles/{}/'.format(group_name)].keys():
        print('Found id/ group, columns will be sorted.')
        ids = h5file['/particles/{}/id/value'.format(group_name)][begin:end:step]

    # Prepares box. Assumes that box is static even if there are time-dependent values.
    box = h5file['/particles/{}/box/edges'.format(group_name)]
    if 'value' in box:
        box = numpy.array(box['value'][0])
    else:
        box = numpy.array(box)

    # Preapres trajectory with image convention.
    trj = numpy.array(
        h5file['/particles/{}/position/value'.format(group_name)][begin:end:step]
        )
    if ids is not None:
        trj = sort_h5md_array(trj, ids)

    if 'image' in h5file['/particles/{}'.format(group_name)].keys() and not no_image:
        print('Found image group, computing absolute trajectory...')
        image = numpy.array(
            h5file['/particles/{}/image/value'.format(group_name)][begin:end:step]
        )
        if ids is not None:
            image = sort_h5md_array(image, ids)
        trj = trj + box*image

    # Prepares masses.
    masses = h5file['/particles/{}/mass'.format(group_name)]
    if 'value' in masses:
        masses = masses['value']
        if ids is not None:
            masses = sort_h5md_array(masses, ids, 1)
        masses = masses[0]
    masses = numpy.array(masses)
    return ids, box, trj, masses

class File(object):
  """File object.

  Args:
    opened: The flag indicated that the file is opened.
    content: The raw content of the file.
    chains: The dict with chains.
    chain_idx_type_map: The map from idx to type.
    chain_type_idx_map: The map from type to chain idx.
    box: The tuple with simulation box dimensions.
    scale_factor: How to scale the numbers, by default we use Angstroms.
  """
  opened = False
  content = None
  chains = {}
  chain_idx_type_map = {}
  chain_type_idx_map = {}
  box = None
  data = None
  scale_factor = 1.0
  file = None


  def __init__(self, file_name, scale_factor=None):
    self.file_name = file_name
    if scale_factor:
        self.scale_factor = scale_factor

  def open(self):
    """Open the pdb file."""
    self.file = open(self.file_name, 'r')
    self.opened = True


class GROFile(File):
  scale_factor = 10.0

  def read(self):
    """Reads the .gro file and return the atom list.

    Returns:
      The dict with atoms (key: atom_id, value: atom object).
    """

    atoms = {}
    if not self.opened:
      raise Exception('File is not open.')

    if not self.content:
      self.content = self.file.readlines()

    number_of_atoms = int(self.content[1])

    for line in self.content[2:number_of_atoms + 2]:
      chain_idx = int(line[0:5].strip())
      chain_name = line[5:10].strip()
      at_name = line[10:15].strip()
      at_id = int(line[15:20].strip())
      # Nedd to rescale.
      pos_x = float(line[20:28].strip()) * self.scale_factor
      pos_y = float(line[28:36].strip()) * self.scale_factor
      pos_z = float(line[36:44].strip()) * self.scale_factor

      atoms[at_id] = (
          structures.Atom(
              atom_id=at_id,
              name=at_name,
              chain_name=chain_name,
              chain_idx=chain_idx,
              position=numpy.array([pos_x, pos_y, pos_z])
          ))

      if chain_idx not in self.chains:
        self.chains[chain_idx] = set()
        self.chain_idx_type_map[chain_idx] = chain_name
      self.chains[chain_idx].add(atoms[at_id])

    # Reads the box size, the last line.
    self.box = numpy.array(
        map(float, filter(None, self.content[number_of_atoms + 2].split(' ')))
        ) * self.scale_factor
    self.data = atoms

  def write(self, file_name=None, force=False):
    """Writes the content to the output file.

    Args:
        file_name: The new file name, otherwise the old one will be used.
        force: Force to save even if any atoms wasn't updated.
    """
    output = []
    output.append('XXX of molecules')
    # Puts the number of atoms
    output.append('%d' % len(self.data))
    # Puts the definition of the atoms, fixed format.
    fmt = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
    for at_id in sorted(self.data):
        at = self.data[at_id]
        output.append(fmt % (
            at.chain_idx,
            at.chain_name,
            at.name,
            at.atom_id,
            at.position[0],
            at.position[1],
            at.position[2]
            ))

    output.append('%f %f %f' % tuple(self.box))
    write_file_path = prepare_path(file_name if file_name else self.file_name)
    print('Writing GRO file %s', write_file_path)
    output_file = open(write_file_path, 'w')
    output_file.writelines('\n'.join(output))


class GROTrajectory(File):
    trajectory = []
    time_steps = []
    box = []
    def read(self):
        if not self.opened:
            raise Exception('File is not open.')
        current_step = 0
        in_header = False
        for line in self.file:
            if not line.strip():
                continue
            if 'current_step' in line or 't=' in line:
                current_step = int(line.split(',')[1].split('=')[1])
                self.time_steps.append(current_step)
                if self.trajectory and self.trajectory[-1]:
                    self.trajectory[-1] = numpy.array(self.trajectory[-1])
                self.trajectory.append([])
                print current_step
                in_header = True
                continue
            elif in_header:
                in_header = False
                continue
            elif len(line.split()) == 3:
                self.box = map(float, line.split())
            else:
                # Nedd to rescale.
                pos_x = float(line[20:28].strip()) * self.scale_factor
                pos_y = float(line[28:36].strip()) * self.scale_factor
                pos_z = float(line[36:44].strip()) * self.scale_factor
                self.trajectory[-1].append([pos_x, pos_y, pos_z])


class PDBFile(File):

  def read(self):
    """Reads the file and return atom list."""
    atoms = {}
    if not self.opened:
      raise Exception('File is not open.')

    if not self.content:
      self.content = self.file.readlines()

    for line in self.content:
      if line.startswith('CRYST1'):
        # Box size
        self.box = numpy.array(
            map(float, filter(None, line.split(' '))[1:4])
            )
      elif line.startswith('ATOM') or line.startswith('HETATM'):
        atom_id = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        chain_name = line[17:20].strip()  # Residue name
        chain_idx = line[22:26].strip()
        pos_x = float(line[30:38])
        pos_y = float(line[38:46])
        pos_z = float(line[46:54])
        atoms[atom_id] = (
          structures.Atom(
              atom_id=atom_id,
              name=atom_name,
              chain_name=chain_name,
              chain_idx=chain_idx,
              position=numpy.array([pos_x, pos_y, pos_z])
              ))
        if chain_idx not in self.chains:
          self.chains[chain_idx] = set()
          self.chain_idx_type_map[chain_idx] = chain_name
        self.chains[chain_idx].add(atoms[atom_id])

    self.data = atoms


  def write(self, output_file, line_replace=None):
    """Write the file again."""
    if not self.opened:
      raise Exception('File is not open.')

    for line in self.content:
      line_replacer = line_replace.get(line[0:6].strip())
      if line_replacer:
        output_file.write(line_replacer(line))
      else:
        output_file.write(line)

  def write_full(self, file_name=None):
    """Write the file again."""
    output = []
    output.append('REMARK generate by YAPT')
    output.append('MODEL 1')
    # Writing the box coordinates
    # Following http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM
    # Boxes are orthorhombic for now
    output.append('%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n' % (
        'CRYST1',
        self.box[0] / self.scale_factor,
        self.box[1] / self.scale_factor,
        self.box[2] / self.scale_factor,
        90.00,
        90.00,
        90,
        'P 1',
        1
        ))

    # Puts the number of atoms
    output.append('%d' % len(self.atoms))
    # Puts the definition of the atoms, fixed format.
    fmt = '%-6s%5d %4s %-3s  %4d    %8.3f%8.3f%8.3f                      %2s'
    for at_id in sorted(self.atoms):
        at = self.atoms[at_id]
        output.append(fmt % (
            'ATOM  ',
            int(at.atom_id) % 100000,
            at.name,
            at.chain_name,
            int(at.chain_idx) % 10000,
            at.position[0] / self.scale_factor,
            at.position[1] / self.scale_factor,
            at.position[2] / self.scale_factor,
            at.name
            ))

    output.append('TER')
    output.append('ENDMDL')

    write_file_path = prepare_path(file_name if file_name else self.file_name)
    output_file = open(write_file_path, 'w')
    output_file.writelines('\n'.join(output))
    output_file.close()
    self.atoms_updated = False


class TopologyFile(File):
  """Very basic representation of the GROMACS topology file."""

  def __init__(self, file_name):
    super(TopologyFile, self).__init__(file_name)
    self.parsers = {
        'atoms': self._parse_atoms,
        'bonds': self._parse_bonds,
        'dihedrals': self._parse_dihedrals,
        'improper_dihedrals': self._parse_improper_dihedrals,
        'pairs': self._parse_pairs,
        'angles': self._parse_angles
        }

    self.writers = {
        'atoms': self._write_atoms,
        'bonds': self._write_bonds,
        'angles': self._write_angles,
        'dihedrals': self._write_dihedrals,
        'improper_dihedrals': self._write_improper_dihedrals,
        'pairs': self._write_pairs
        }
    # Current set of data.
    self.angles = {}
    self.atoms = {}
    self.bonds = {}
    self.chains = {}
    self.chains_atoms = {}
    self.chains_idx = {}
    self.chain_idx = {}
    self.dihedrals = {}
    self.improper_dihedrals = {}
    self.pairs = {}

    # Store the new data.
    self.new_data = {
        'atoms': {},
        'bonds': {},
        'angles': {},
        'dihedrals': {},
        'improper_dihedrals': {},
        'pairs': {}
        }
    # Some helper structures for bonds.
    # key: atom_i
    # val: set of atom_j
    self.bonds_def = defaultdict(set)


  def update_position(self, pdbfile):
    """Reads the position data from the coordinate file and update the atoms.

    Args:
      pdbfile: The pdb file.
    """
    for k, v in pdbfile.data.iteritems():
      self.atoms[k].position = v.position


  def add_bonds(self, new_bonds):
    """Adds the new bonds to the current set.

    Args:
      new_bonds: The dict with the key as atom_pair and value as the list
        of new bonds
    """
    for bonds in new_bonds.values():
      for bond in bonds:
        self.bonds.append(bond)
        if bond[0] not in self.bonds_def:
          self.bonds_def[bond[0]] = set()
        self.bonds_def[bond[0]].add(bond[1])
        if bond[1] not in self.bonds_def:
          self.bonds_def[bond[1]] = set()
        self.bonds_def[bond[1]].add(bond[0])

  def remove_bonds(self, bond_list):
    """Remove bonds and all linked information within.

    Args:
      bond_list: The list of tuples with bond definitions.
    """

    for b1, b2 in bond_list:
      # remove bond definition
      b_set = set([b1, b2])
      bonds_def = False
      if (b1, b2) in self.bonds:
        del self.bonds[(b1, b2)]
        bonds_def = True
      elif (b2, b1) in self.bonds:
        del self.bonds[(b2, b1)]
        bonds_def = True

      if (b1, b2) in self.new_data['bonds']:
        del self.new_data['bonds'][(b1, b2)]
        bonds_def = True
      elif  (b2, b1) in self.new_data['bonds']:
        del self.new_data['bonds'][(b2, b1)]
        bonds_def = True

      if bonds_def:
        self.bonds_def[b1].remove(b2)
        self.bonds_def[b2].remove(b1)

      # Clean angles
      for ang_def in self.angles:
        if not b_set - set(ang_def):
          del self.angles[ang_def]
      for ang_def in self.new_data['angles']:
        if not b_set - set(ang_def):
          del self.new_data['angles'][ang_def]

      # Clean dihedrals
      for dih_def in self.dihedrals:
        if not b_set - set(dih_def):
          del self.dihedrals[dih_def]
      for dih_def in self.new_data['dihedrals']:
        if not b_set - set(dih_def):
          del self.new_data['dihedrals'][dih_def]

      for dih_def in self.improper_dihedrals:
        if not b_set - set(dih_def):
          del self.improper_dihedrals[dih_def]
      for dih_def in self.new_data['improper_dihedrals']:
        if not b_set - set(dih_def):
          del self.new_data['improper_dihedrals'][dih_def]

      # Clean pairs
      if (b1, b2) in self.pairs:
        del self.pairs[(b1, b2)]
      elif (b2, b1) in self.pairs:
        del self.pairs[(b2, b1)]

      if (b1, b2) in self.new_data['pairs']:
        del self.new_data['pairs'][(b1, b2)]
      elif (b2, b1) in self.new_data['pairs']:
        del self.new_data['pairs'][(b2, b1)]

  # Parsers for the data.
  def _parse_bonds(self, raw_data):
    atom_tuple = tuple(map(int, raw_data[0:2]))
    self.bonds[atom_tuple] = raw_data[2:]

    self.bonds_def[atom_tuple[0]].add(atom_tuple[1])
    self.bonds_def[atom_tuple[1]].add(atom_tuple[0])

  def _parse_atoms(self, raw_data):
    at = structures.Atom(
        atom_id=int(raw_data[0]),
        atom_type=raw_data[1],
        chain_idx=int(raw_data[2]),
        chain_name=raw_data[3],
        name=raw_data[4],
        cgnr=int(raw_data[5])
        )

    if len(raw_data) > 6:
      at.charge = float(raw_data[6])
    if len(raw_data) > 7:
      at.mass = float(raw_data[7])

    if at.chain_name not in self.chains:
      self.chains[at.chain_name] = defaultdict(list)
      self.chains_atoms[at.chain_name] = defaultdict(dict)

    if at.chain_idx not in self.chains_idx:
      self.chains_idx[at.chain_idx] = set()
    self.chains_idx[at.chain_idx].add(at)

    if at.chain_name not in self.chain_idx:
      self.chain_idx[at.chain_name] = {}
    if at.chain_idx not in self.chain_idx[at.chain_name]:
      self.chain_idx[at.chain_name][at.chain_idx] = {}

    self.chains[at.chain_name][at.chain_idx].append(at)
    self.chain_idx[at.chain_name][at.chain_idx][at.atom_id] = at

    self.atoms[at.atom_id] = at

  def _parse_angles(self, raw_data):
    atom_tuple = tuple(map(int, raw_data[0:3]))
    self.angles[atom_tuple] = raw_data[3:]

  def _parse_dihedrals(self, raw_data):
    atom_tuple = tuple(map(int, raw_data[0:4]))
    self.dihedrals[atom_tuple] = raw_data[4:]

  def _parse_improper_dihedrals(self, raw_data):
    atom_tuple = tuple(map(int, raw_data[0:4]))
    self.improper_dihedrals[atom_tuple] = raw_data[4:]

  def _parse_pairs(self, raw_data):
    atom_tuple = tuple(map(int, raw_data[0:2]))
    self.pairs[atom_tuple] = raw_data[2:]

  # Writers
  def _write_atoms(self):
    return_data = []
    def __write_atoms(atoms_dict):
      r = []
      for at_id in sorted(atoms_dict):
        at = atoms_dict[at_id]
        at_line = ('%s ' * len(at)).strip()
        r.append(at_line % tuple(at))
      return r

    if self.atoms:
      atoms_dict = {k: v.get_tuple() for k, v in self.atoms.iteritems()}
      return_data.extend(__write_atoms(atoms_dict))
    if self.new_data['atoms']:
      return_data.append('; new atoms')
      return_data.extend(__write_atoms(self.new_data['atoms']))

    return return_data

  def _write_bonds(self):  # pylint:disable=R0201
    return_data = []
    if self.bonds:
      return_data.extend(self._write_default(self.bonds))
    if self.new_data['bonds']:
      return_data.append('; new bonds')
      for key, values in self.new_data['bonds'].iteritems():
        if (key[0], key[1]) not in self.bonds or (key[1], key[0]) not in self.bonds:
          return_data.append('%d  %d  %s' % (
            key[0], key[1], ' '.join(map(str, values))
            ))
    return return_data

  def _write_pairs(self):  # pylint:disable=R0201
    return_data = []
    if self.pairs:
      return_data.extend(self._write_default(self.pairs))
    if self.new_data['pairs']:
      return_data.append('; new pairs')
      return_data.extend(self._write_default(self.new_data['pairs']))
    return return_data

  def _write_angles(self):
    return_data = []
    if self.angles:
      return_data.extend(self._write_default(self.angles))
    if self.new_data['angles']:
      return_data.append('; new angles')
      return_data.extend(self._write_default(self.new_data['angles']))
    return return_data

  def _write_dihedrals(self):
    return_data = []
    if self.dihedrals:
      return_data.extend(self._write_default(self.dihedrals))
    if self.new_data['dihedrals']:
      return_data.append('; new dihedrals')
      return_data.extend(self._write_default(self.new_data['dihedrals']))
    return return_data

  def _write_improper_dihedrals(self):
    return_data = []
    if self.improper_dihedrals:
      return_data.extend(self._write_default(self.improper_dihedrals))
    if self.new_data['improper_dihedrals']:
      return_data.append('; new improper dihedrals')
      return_data.extend(self._write_default(self.new_data['improper_dihedrals']))
    return return_data


  def _write_default(self, data):  # pylint:disable=R0201
    flat_data = [list(k) + list(v) for k, v in data.iteritems()]
    flat_data.sort()
    return ['%s' % ' '.join(map(str, x)) for x in flat_data]

  def read(self):
    """Reads the topology file."""

    self.file = open(self.file_name, 'r')

    if not self.content:
      self.content = self.file.readlines()

    # New version
    current_parser = None
    visited_sections = set()
    section_name = None
    previous_section = None
    for line in self.content:
      line = line.strip()
      if line.startswith(';') or line.startswith('#') or len(line) == 0:
        continue
      elif line.startswith('['):  # Section
        previous_section = section_name
        section_name = line.replace('[', '').replace(']', '').strip()
        # Hack for GROMACS improper_dihedrals
        if previous_section == 'dihedrals' and section_name == 'dihedrals':
          section_name = 'improper_dihedrals'
        current_parser = self.parsers.get(section_name)
        visited_sections.add(previous_section)
      else:
        if current_parser is not None and section_name not in visited_sections:
          raw_data = filter(None, line.split())
          if raw_data:
            current_parser(raw_data)  # pylint:disable=E1102

    self._postread()

  def _postread(self):
    """Runs after the reading sections."""
    # Updates the degree.
    for b1, b2 in self.bonds:
      self.atoms[b1].degree += 1
      self.atoms[b2].degree += 1

    # Updates the neighbours.
    for at in self.atoms.values():
      at.neighbours = [
          self.atoms[j] for j in self.bonds_def.get(at.atom_id, [])
          ]

  def write(self, filename=None):
    """Updates the topology file.

    Args:
      filename: The optional output filename.
    """

    if filename:
      output_file = open(prepare_path(filename), 'w')
    else:
      output_file = open(prepare_path(self.file_name), 'w')

    new_data = []
    previous_section = None
    current_section = None
    skip_lines = False
    section_writer = None

    if not self.content:
      self.content = [
          '[ atoms ]\n',
          '\n',
          '[ bonds ]\n',
          '\n',
          '[ angles ]\n',
          '\n',
          '[ dihedrals ]\n',
          '\n',
          '[ dihedrals ]\n',
          '\n',
          '[ pairs ]\n',
          '\n'
          ]

    for line in self.content:
      tmp_line = line.strip()
      if tmp_line.startswith('['):  # section part
        new_data.append(line)
        previous_section = current_section
        current_section = tmp_line.replace('[', '').replace(']', '').strip()
        if previous_section == 'dihedrals' and current_section == 'dihedrals':
          current_section = 'improper_dihedrals'
        section_writer = self.writers.get(current_section)
        skip_lines = False
      else:
        if section_writer is None:  # there is no special writer, simply copy the line
          new_data.append(line)
        elif skip_lines == False:
          new_data.extend(['%s\n' % x for x in section_writer()])  # pylint:disable=E1102
          skip_lines = True

    output_file.writelines(new_data)
    output_file.close()


def readtrj(filename):
    return {'gro': GROFile, 'pdb': PDBFile}[filename.split('.')[-1]](filename, scale_factor=1.0)
