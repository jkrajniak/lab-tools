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

BondMatcher - data structures."""

import inspect
import sys

class Atom(object):
  """Define the atom."""
  atom_id = None
  name = None
  chain_name = None
  chain_ids = None
  position = None
  _degree = 0
  max_degree = None
  neighbours = None
  radius = None
  cgnr = None
  charge = None
  mass = None

  def __init__(self,
      atom_id=None,
      atom_type=None,
      name=None,
      chain_name=None,
      position=None,
      degree=0,
      max_degree=None,
      chain_idx=None,
      radius=None,
      charge=None,
      cgnr=None,
      mass=None):
    self.atom_id = atom_id
    self.atom_type = atom_type
    self.name = name
    self.chain_name = chain_name  # Chain name or other id
    self.position = position
    self.radius = radius
    self._degree = degree
    self.max_degree = max_degree
    self.chain_idx = chain_idx  # Sequence number of the chain
    self.charge = charge
    self.cgnr = cgnr
    self.mass = mass

  @staticmethod
  def create(base_atom, atom_id, chain_id):
    return Atom(
        atom_id=atom_id,
        name=base_atom.name,
        chain_name=base_atom.chain_name,
        chain_idx=chain_id,
        max_degree=base_atom.max_degree,
        radius=base_atom.radius,
        charge=base_atom.charge)

  @property
  def degree(self):
    return self._degree

  @degree.setter
  def degree(self, value):
    if self.max_degree is not None and value > self.max_degree:
      raise ValueError('Invalid degree, %d is greater than %d' % (value, self.max_degree))
    self._degree = value

  def __hash__(self):
    # BUG: watch out on the set() with those objects...
    return hash((self.name, self.chain_name))

  def __eq__(self, other):
    """Compares the name. Then you can simple do sth like that:
    'H1' in [AtomDef('H1'), AtomDef('N1')]
    """
    if other is None:
      return False
    return self.chain_name == other.chain_name and self.name == other.name

  def cmp_exact(self, other):
    """Compares chain name, name, chain id, atom id."""
    if other is None:
      return False

    return (self.chain_name == other.chain_name and self.name == other.name
            and self.atom_id == other.atom_id and self.chain_idx == other.chain_idx
           )

  def __repr__(self):
    return '(id=%s, name=%s, chain_name=%s, position=%s, degree=%s, radius=%s, charge=%s)' % (
        self.atom_id, self.name, self.chain_name, self.position, self.degree,
        self.radius, self.charge
        )

  def __unicode__(self):
    return 'A%s' % self.__repr__()

  def get_dict(self):
    """Returns the dictionary."""
    return {
        'id': self.atom_id,
        'name': self.name,
        'chain_name': self.chain_name,
        'chain_idx': self.chain_idx,
        'position': self.position,
        'degree': self._degree,
        'max_degree': self.max_degree,
        'radius': self.radius,
        'charge': self.charge
        }

  @property
  def sstr(self):
    """Returns simple string."""
    return (self.name, self.chain_name)

  def get_tuple(self):
    """Returns tuple."""
    ret = [
        self.atom_id,
        self.atom_type,
        self.chain_idx,
        self.chain_name,
        self.name,
        self.cgnr
        ]
    if self.charge:
      ret.append(self.charge)
    if self.mass:
      if not self.charge:
        ret.append(0.0)
      ret.append(self.mass)

    return ret


class BaseSettings(object):
  def __init__(self, settings_file=None):
    settings_variables = {}
    if settings_file:
      copy_globals = globals().copy()
      try:
        execfile(settings_file, copy_globals, settings_variables)
      except SyntaxError as ex:
        print 'There is an error in your config file %s in line %s' % (ex.filename, ex.lineno)
        print '%d: %s' % (ex.lineno, ex.text)
        print ex.msg
        sys.exit(2)

      set_settings_variables = set(settings_variables)
      missing_sections = self._required_sections - set_settings_variables

      if missing_sections:
        raise Exception(
          'Invalid settings file, those sections are missing: %s' % ','.join(
            missing_sections))

      obsolete_sections = self._obsolete_sections.intersection(set_settings_variables)
      if obsolete_sections:
        print 'Those sections are obsolete: %s' % ','.join(obsolete_sections)

      replace_sections = set(self._replace_sections).intersection(set_settings_variables)
      if replace_sections:
        print 'Warning! Please use those section names:\n%s' % (
            '\n'.join(['- %s instead of %s' % (x, self._replace_sections[x])
                       for x in replace_sections])
            )

        # Replace the section names.
        for k, v in self._replace_sections.iteritems():
          if k in settings_variables:
            settings_variables[v] = settings_variables[k]

      # Load the variables as it they would be and attributes of the class.
      self.__dict__.update(settings_variables)


# Base class for settings.
class BondSettings(BaseSettings):
  _important_atoms = None
  INTRAMOLECULAR_BONDS = False
  NEIGHBOURS = {}

  _required_sections = set([
      'ANGLES',
      'DIHEDRALS',
      'PAIRS',
      'ATOM_PAIRS',
      'BOND_CONFIG',
      'INIT_CUTOFF_DISTANCE',
      'CUTOFF_STEP'
      ])
  _obsolete_sections = set([])
  _replace_sections = {}

  def __init__(self, settings_file=None):
    super(BondSettings, self).__init__(settings_file)

    # Search for the important atoms.
    self.IMPORTANT_ATOMS = [
        x[1] for x in inspect.getmembers(self) if isinstance(x[1], Atom)
        ]
    if not self.IMPORTANT_ATOMS:
      raise Exception('Atom definition not found.')

    self.ATOMS_DEF = {
        (x.name, x.chain_name): x for x in self.IMPORTANT_ATOMS
        }
    self.process_settings()

  def process_settings(self):
    """Generate settings entries that based on existing once."""

    angles_nb_index = {
        0: [1],
        1: [0, 2]
        }
    dihedrals_nb_index = {
        0: [1],
        1: [0, 2],
        2: [1, 3]
        }
    # Lets create the neighbors config based on the ANGLES, DIHEDRALS, PAIRS setting
    for angle_def in self.ANGLES.values():
      for atom_set in angle_def:
        for ai, aj_set in angles_nb_index.iteritems():
          if atom_set[ai] not in self.NEIGHBOURS:
            self.NEIGHBOURS[atom_set[ai]] = set()
          for aj in aj_set:
            self.NEIGHBOURS[atom_set[ai]].add(atom_set[aj])
    for angle_def in self.DIHEDRALS.values():
      for atom_set in angle_def:
        for ai, aj_set in dihedrals_nb_index.iteritems():
          if atom_set[ai] not in self.NEIGHBOURS:
            self.NEIGHBOURS[atom_set[ai]] = set()
          for aj in aj_set:
            self.NEIGHBOURS[atom_set[ai]].add(atom_set[aj])
