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

from md_libs import files_io


def _args():
    parser = argparse.ArgumentParser('Update charges in GROMACS topology')
    parser.add_argument('input_topol', help='Input topology')
    parser.add_argument('charge_map', help='New charges, format: atom_name:charge,atom_name:charge')
    parser.add_argument('--zero_charge', action='store_true',
                        help=('If set to true, then the total charge of molecules'
                              ' in topology is made to be zero'))

    return parser.parse_args()


def main():
    args = _args()

    in_top = files_io.GROMACSTopologyFile(args.input_topol)
    in_top.read()

    at_name2charge = {}
    for x in args.charge_map.split(','):
        z = x.split(':')
        at_name2charge[z[0]] = float(z[1])

    print(at_name2charge)
    current_charge = sum([x.charge for x in in_top.atoms.values()])
    print('Current total charge: {}'.format(current_charge))

    for at_data in in_top.atoms.values():
        if at_data.name in at_name2charge:
            at_data.charge = at_name2charge[at_data.name]
    current_charge = sum([x.charge for x in in_top.atoms.values()])
    print('New total charge: {}'.format(current_charge))
    if args.zero_charge:
        charge_map = {x.name: x.charge for x in in_top.atoms.values()}
        num_charges = len(in_top.atoms) - len([k for k in at_name2charge if k in charge_map])
        delta_charge = current_charge / num_charges
        for at_data in in_top.atoms.values():
            if at_data.name not in at_name2charge:
                at_data.charge -= delta_charge
        current_charge = sum([x.charge for x in in_top.atoms.values()])
        print('New total charge: {}'.format(current_charge))

    in_top.write(force=True)

if __name__ == '__main__':
    main()

