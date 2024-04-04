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
import h5py


def _args():
    parser = argparse.ArgumentParser("Basic information about h5md trajectory file.")
    parser.add_argument("input_file")

    return parser.parse_args()


def main():
    args = _args()

    h5 = h5py.File(args.input_file, "r")

    groups = list(h5["/particles/"].keys())
    max_time = max(h5["/particles/{}/position/time".format(groups[0])])
    print(("File: {}".format(args.input_file)))
    print(("H5MD groups: {}".format(groups)))
    print(("Max time: {}".format(max_time)))

    h5.close()


if __name__ == "__main__":
    main()
