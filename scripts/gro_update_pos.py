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
    parser = argparse.ArgumentParser("Copy position from one gro file to other")
    parser.add_argument("input_gro")
    parser.add_argument("positions")
    parser.add_argument("output_gro")

    return parser.parse_args()


def main():
    args = _args()

    input_gro = files_io.GROFile(args.input_gro)
    input_gro.read()
    positions_gro = files_io.GROFile(args.positions)
    positions_gro.read()

    # Update positions, based on atom id.
    print("Update positions, based on atom_id")
    for at_id, at_data in list(positions_gro.atoms.items()):
        input_gro.atoms[at_id] = input_gro.atoms[at_id]._replace(position=at_data.position)

    input_gro.write(args.output_gro, force=True)


if __name__ == "__main__":
    main()
