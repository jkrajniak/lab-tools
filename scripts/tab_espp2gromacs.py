#!/usr/bin/env python
"""
Copyright (C) 2015 Jakub Krajniak <jkrajniak@gmail.com>

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
import datetime
import numpy as np


def _args():
    parser = argparse.ArgumentParser(description="Converts Espresso++ to GROMACS")
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    parser.add_argument("--table_type", choices=("pair", "bond", "angle", "dihedral"), default="pair")
    parser.add_argument("--length_scale", default=1.0, type=float)

    return parser.parse_args()


def _pair_convert(input_f, output_f, args):
    in_f = np.loadtxt(input_f)
    out_f = np.zeros((in_f.shape[0], 7))
    out_f[:, 0] = in_f[:, 0]
    out_f[:, 5] = in_f[:, 1]
    np.savetxt(
        output_f, out_f, header="Converted ESPP to GROMACS, {} to {}, {}".format(input_f, output_f, datetime.datetime.now())
    )


def main():
    args = _args()

    table_type2func = {
        "pair": _pair_convert,
    }

    table_type2func[args.table_type](args.input_file, args.output_file, args)


if __name__ == "__main__":
    main()
