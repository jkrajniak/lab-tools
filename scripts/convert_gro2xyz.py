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

from md_libs import files_io

import argparse
import re


def _args():
    parser = argparse.ArgumentParser("Convert .gro file to XYZ file format")
    parser.add_argument("gro", help="Input .gro file")
    parser.add_argument("xyz", help="Output .xyz file")
    parser.add_argument("--comment_line", default="")
    parser.add_argument("--name_map", help="Convert atom names to atom symbols, e.g. C1->C; format: C1:C,C2:C,C3:C,H11:H")
    parser.add_argument(
        "--use_first_alphabetic",
        help="Convert atom names to atom symbols by using only letters, e.g. C1->C, H11->H",
        action="store_true",
    )

    return parser.parse_args()


def main():
    args = _args()

    if args.name_map and args.use_first_alphabetic:
        print("You can use either --name_map or --use_first_alphabetic but not both")
        return False

    name_map = {}
    if args.name_map:
        name_map = dict([k.split(":") for k in args.name_map.split(",")])

    only_alpha = re.compile("[^a-zA-Z]+")

    gro_file = files_io.GROFile(args.gro)
    gro_file.read()

    xyz_file = files_io.XYZFile(args.xyz)
    for at_id in gro_file.atoms:
        pos = gro_file.atoms[at_id].position
        at_name = gro_file.atoms[at_id].name
        if args.name_map:
            at_name = name_map.get(at_name, at_name)
        elif args.use_first_alphabetic:
            at_name = re.sub(only_alpha, "", at_name)
        xyz_file.atoms[at_id] = gro_file.atoms[at_id]._replace(position=10.0 * pos, name=at_name)

    xyz_file.write()


if __name__ == "__main__":
    main()
