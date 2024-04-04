#!/usr/bin/env python
"""
Copyright (C) 2019 Jakub Krajniak <jkrajniak@gmail.com>

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
    parser = argparse.ArgumentParser("Synchronize the atom and chain names from topology to gro file")
    parser.add_argument("--input_gro", required=True)
    parser.add_argument("--input_top", required=True)
    parser.add_argument("--output_gro", required=True)
    return parser.parse_args()


def main():
    args = _args()

    in_gro = files_io.GROFile(args.input_gro)
    in_top = files_io.GROMACSTopologyFile(args.input_top)
    in_gro.read()
    in_top.read()

    import IPython

    IPython.embed()

    for a_id, a in list(in_top.atoms.items()):
        in_gro.atoms[a_id] = in_gro.atoms[a_id]._replace(chain_name=a.chain_name, name=a.name)

    in_gro.write(args.output_gro, force=True)


if __name__ == "__main__":
    main()
