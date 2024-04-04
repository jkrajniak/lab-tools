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
    parser = argparse.ArgumentParser("Read GAUSSIAN log file and copy partial charges to GROMASC topology")
    parser.add_argument("gaussian_log", help="Input GAUSSIAN log file")
    parser.add_argument("input_topology", help="Input GROMACS topology")

    return parser.parse_args()


def main():
    args = _args()
    in_esp_charges = False
    espl = []
    with open(args.gaussian_log, "r") as gaussian_log:
        for l in gaussian_log:
            if "ESP charges:" in l:
                in_esp_charges = True
                espl = []
            elif "Sum of ESP charges" in l:
                in_esp_charges = False
            elif in_esp_charges:
                tmpl = l.split()
                print(tmpl)
                if len(tmpl) > 1:
                    espl.append((tmpl[1], tmpl[2]))

    print(espl)

    in_top = files_io.GROMACSTopologyFile(args.input_topology)
    in_top.read()

    for at_idx, at_id in enumerate(sorted(in_top.atoms)):
        print(at_idx, espl[at_idx], in_top.atoms[at_id].name)
        at_prefix, new_charge = espl[at_idx]
        if in_top.atoms[at_id].name.startswith(at_prefix):
            if in_top.atoms[at_id].charge != float(new_charge):
                print(("New charge {} old charge {}".format(new_charge, in_top.atoms[at_id].charge)))
            in_top.atoms[at_id].charge = float(new_charge)
        else:
            raise RuntimeError("Wrong atom prefix")
    in_top.write(force=True)


if __name__ == "__main__":
    main()
