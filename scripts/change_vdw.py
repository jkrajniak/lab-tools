#!/usr/bin/env python
"""
Copyright (C) 2014 Jakub Krajniak <jkrajniak@gmail.com>

This program is free software: you can redistribute it and/or modify
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


parser = argparse.ArgumentParser(description="Change the vdw radius")
parser.add_argument("-i", "--input_itp", help="Input itp file", required=True)
parser.add_argument("-o", "--output_itp", help="Output itp file", required=True)
parser.add_argument("-c", "--config", help="Config file", required=True)
parser.add_argument("-p", "--percentage", help="The VdW radius in percanteg [0.0-1.0]", type=float)

args = parser.parse_args()

exec(compile(open(args.config, "rb").read(), args.config, "exec"))

input_data = open(args.input_itp, "r").readlines()
output_file = open(args.output_itp, "w")

if args.percentage < 0.0 or args.percentage > 1.0:
    raise ValueError("The percentage has to be in the range 0.0-1.0")

atom_types = {}
section = None
new_data = []
for line in input_data:
    if line.startswith("["):
        section = line.strip().replace("[", "").replace("]", "").strip()
        if section != "atomtypes":
            new_data.append(line)
    else:
        if section == "atomtypes":
            tt = list(map(str.strip, line.split()))
            if len(tt) > 0:
                atom_types[tt[0]] = tt
        else:
            new_data.append(line)

for atom_type, atom_def in atom_types.items():
    # Let's change the VdW radius.
    atom_def[6] = "%e" % (args.percentage * DEFAULT_VDW[atom_type][0])

# Write the new file
output_file.write("[ atomtypes ]\n")
output_file.writelines(["%s\n" % "\t".join(x) for x in list(atom_types.values())])
output_file.write("\n")
output_file.writelines(new_data)
output_file.close()
