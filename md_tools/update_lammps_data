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

__doc__ = 'Update lammps data file from other data file'


def _args():
    parser = argparse.ArgumentParser(
        description='Update LAMMPS data file',
        add_help=True)

    parser.add_argument('--in_data', help='Inpput data file', required=True)
    parser.add_argument('--out_data', help='Output LAMMPS data file', required=True)
    parser.add_argument('--what_to_update', help='What to update, colon separate list: coeff,position')

    return parser

def main():
    args = _args().parse_args()

    input_data = files_io.LammpsReader()
    input_data.read_data(args.in_data)




if __name__ == '__main__':
    main()