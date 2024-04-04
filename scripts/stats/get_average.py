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
import numpy as np


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("column", type=int)
    parser.add_argument("--rescale", type=float, default=1.0)
    parser.add_argument("--block_size", type=float, default=22.0)
    parser.add_argument("--append", default=None)

    return parser.parse_args()


def main():
    args = _args()
    input_data = np.loadtxt(args.input_file)[:, args.column] * args.rescale

    avg = np.average(input_data[:: int(args.block_size)])
    std = np.std(input_data) * np.sqrt(args.block_size / input_data.shape[0])

    if args.append:
        with open(args.append, "a+") as outf:
            outf.write("{} {}\n".format(avg, std))
    else:
        print("avg: {} std:{}".format(avg, std))


if __name__ == "__main__":
    main()
