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
from md_libs import files_io


def _args():
    parser = argparse.ArgumentParser("Check if the atom labels are unique.")
    parser.add_argument("--coord", help="Input .gro coordinate file")
    parser.add_argument("--top", help="Input .top file")

    return parser.parse_args()


def check_file(input_file, file_name):
    atom_symbols = [x.name for x in list(input_file.atoms.values())]
    atom_sym_count = [(x, atom_symbols.count(x)) for x in atom_symbols]
    file_correct = True
    for sym, cn in atom_sym_count:
        if cn > 1:
            file_correct = False
            print(("{}: Symbol {} found {} times".format(file_name, sym, cn)))

    return file_correct


def main():
    args = _args()

    if args.coord:
        input_file = files_io.GROFile(args.coord)
        input_file.read()
        if check_file(input_file, args.coord):
            print(("File {} is correct".format(args.coord)))

    if args.top:
        input_file = files_io.GROMACSTopologyFile(args.top)
        input_file.read()
        total_charge = sum([x.charge for x in list(input_file.atoms.values())])
        if check_file(input_file, args.top):
            print(("File {} is correct".format(args.top)))
        print(("{}: total charge {}".format(args.top, total_charge)))


if __name__ == "__main__":
    main()
