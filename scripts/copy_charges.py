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
    parser = argparse.ArgumentParser("Copy charges from source GROMASC topology to output")
    parser.add_argument("source_topol", help="Source topology file")
    parser.add_argument("output_topol", help="Output topology file")
    parser.add_argument("--ignore_missing", action="store_true", help="Ignore if atom is missing in source topology")
    parser.add_argument(
        "--value_missing",
        default=None,
        type=float,
        required=False,
        help="Set value for atom that is missing in source topology",
    )

    return parser.parse_args()


def main():
    args = _args()

    input_top = files_io.GROMACSTopologyFile(args.source_topol)
    input_top.read()

    output_top = files_io.GROMACSTopologyFile(args.output_topol)
    output_top.read()

    source_at_name2at_data = {at_data.name: at_data for at_data in list(input_top.atoms.values())}
    output_at_name2at_data = {at_data.name: at_data for at_data in list(output_top.atoms.values())}

    for at_name in output_at_name2at_data:
        source_data = source_at_name2at_data.get(at_name)
        if source_data is None and not args.ignore_missing:
            raise RuntimeError("Missing atom name {} in source topology".format(at_name))
        if source_data is not None:
            output_at_name2at_data[at_name].charge = source_data.charge
        elif args.value_missing is not None:
            output_at_name2at_data[at_name].charge = args.value_missing
            print(("Set value {} for missing atom {}".format(args.value_missing, at_name)))

    output_top.write(force=True)


if __name__ == "__main__":
    main()
