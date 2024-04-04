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
import math

from md_libs import files_io


def _args():
    parser = argparse.ArgumentParser("Creates AdResS index particles")
    parser.add_argument("--conf", required=True)
    parser.add_argument("--hy", type=float, required=True)
    parser.add_argument("--ex", type=float, required=True)
    parser.add_argument("--center", required=True, type=float)
    parser.add_argument("--index_translate", default=0, type=int)

    parser.add_argument("--output", required=True)

    return parser.parse_args()


def write(output_file, data, header):
    n_row = 20
    output_file.write("[ %s ]\n" % header)
    rows = int(math.ceil(len(data) / n_row)) + 1
    for n in range(rows):
        row = data[n * n_row : n * n_row + n_row]
        output_file.write(" ".join(map(str, row)))
        output_file.write("\n")


def main():
    args = _args()

    conf = files_io.readtrj(args.conf)
    conf.open()
    conf.read()

    adress_center = args.center
    at_size = args.ex
    hy_size = args.hy
    box_size = conf.box[0]
    at_range = (adress_center - at_size, adress_center + at_size)
    hy_range = (
        (adress_center - at_size - hy_size, adress_center - at_size),
        (adress_center + at_size, adress_center + at_size + hy_size),
    )
    cg_range = ((0.0, adress_center - at_size - hy_size), (adress_center + at_size + hy_size, box_size))

    print("box_size", box_size, "at_size", at_size, "hy_size", hy_size)
    print("at_range", at_range, "hy_range", hy_range, "cg_range", cg_range)

    pos_at = []
    pos_hy = []
    pos_hy_1 = []
    pos_hy_2 = []
    pos_cg = []
    pos_cg_1 = []
    pos_cg_2 = []
    pos_at_cg = []
    pos_at_cg_1 = []
    pos_at_cg_2 = []

    # half at and half cg
    scheme_1 = []
    scheme_1_range = (
        (adress_center + at_size / 2.0, adress_center + at_size),
        (adress_center + at_size + hy_size, adress_center + at_size + hy_size + at_size / 2.0),
    )

    for at_id, atom in conf.data.items():
        pos = atom.position[0]
        x = atom.atom_id + args.index_translate
        # Special cases
        if (pos >= scheme_1_range[0][0] and pos <= scheme_1_range[0][1]) or (
            pos >= scheme_1_range[1][0] and pos <= scheme_1_range[1][1]
        ):
            scheme_1.append(x)

        if pos >= at_range[0] and pos <= at_range[1]:
            pos_at.append(x)
            pos_at_cg.append(x)
            if pos < adress_center:
                pos_at_cg_1.append(x)
            elif pos > adress_center:
                pos_at_cg_2.append(x)
        elif (pos >= hy_range[0][0] and pos < hy_range[0][1]) or (pos > hy_range[1][0] and pos <= hy_range[1][1]):
            pos_hy.append(x)
            if pos >= hy_range[0][0] and pos < hy_range[0][1]:
                pos_hy_1.append(x)
            elif pos > hy_range[1][0] and pos <= hy_range[1][1]:
                pos_hy_2.append(x)
        else:
            pos_cg.append(x)
            if pos > 0.0 and pos < cg_range[0][1]:
                pos_cg_1.append(x)
                pos_at_cg_1.append(x)
            elif pos > cg_range[1][0]:
                pos_cg_2.append(x)
                pos_at_cg_2.append(x)

    assert (len(pos_at) + len(pos_hy) + len(pos_cg)) == len(conf.data)

    output = open(args.output, "w")

    write(output, pos_at, "AT")
    write(output, pos_hy, "HY")
    write(output, pos_hy_1, "HY_1")
    write(output, pos_hy_2, "HY_2")
    write(output, pos_cg, "CG")
    write(output, pos_cg_1, "CG_1")
    write(output, pos_cg_2, "CG_2")
    write(output, scheme_1, "SCHEME_1")

    print(("Wrote file {}".format(args.output)))
    output.close()


if __name__ == "__main__":
    main()
