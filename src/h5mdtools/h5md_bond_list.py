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
import h5py
import numpy as np
import pickle


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5", help="H5MD file")
    parser.add_argument("--time_frame", help="Time frame", default=-1, type=int)
    parser.add_argument("--out", help="Output pickle", required=True)

    return parser.parse_args()


def main():
    args = _args()

    h5 = h5py.File(args.h5, "r")
    if "connectivity" in h5:
        output_connectivities = {}
        for name, ds in list(h5["/connectivity"].items()):
            if isinstance(ds, h5py.Dataset):
                print(("Reading {}".format(name)))
                output_connectivities[name] = np.array(ds)
            else:
                print(
                    (
                        "Reading {}, time frame: {} of {}".format(
                            name, "last" if args.time_frame == -1 else args.time_frame, ds["step"].shape[0]
                        )
                    )
                )
                data = ds["value"][args.time_frame]
                output_connectivities[name] = np.array(data)
        if output_connectivities:
            print(("Writing data to {}".format(args.out)))
            out_file = open(args.out, "wb")
            pickle.dump(output_connectivities, out_file)
            out_file.close()


if __name__ == "__main__":
    main()
