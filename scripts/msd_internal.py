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
import functools
from libs import bonds
import h5py
import multiprocessing as mp
import numpy


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nchains", required=True, type=int)
    parser.add_argument("--chain_length", required=True, type=int)
    parser.add_argument("--begin", default=0, type=int)
    parser.add_argument("--end", default=-1, type=int)
    parser.add_argument("--output")
    parser.add_argument("--output_prefix", default="")
    parser.add_argument("--csv", action="store_true")
    parser.add_argument("--group", default="atoms")
    parser.add_argument("--nt", default=1, type=int)
    parser.add_argument("in_file")

    return parser.parse_args()


def calculate_msd_int(nchains, chain_length, box, trj):
    return bonds.calculate_msd_internal_distance(trj, nchains, chain_length, box, 0.5 * box)


if __name__ == "__main__":
    args = _args()
    data = h5py.File(args.in_file)

    trj = data["/particles/{}/position/value".format(args.group)][args.begin : args.end]
    box = numpy.array(data["/particles/{}/box/edges".format(args.group)])
    if "image" in list(data["/particles/{}".format(args.group)].keys()):
        print("Found image dataset, computing absolute trajectory...")
        image = numpy.array(data["/particles/{}/image/value".format(args.group)][args.begin : args.end])
        trj = trj + box * image

    nt = args.nt
    args.begin = 0
    if args.end == -1:
        args.end = len(trj)
    num_frames = args.end
    print(("Trajectory length: {}".format(num_frames)))
    print(("Number of chains: {}".format(args.nchains)))
    print(("Length of chain: {}".format(args.chain_length)))

    print(("Distribute work on {} process".format(args.nt)))
    func = functools.partial(calculate_msd_int, args.nchains, args.chain_length, box)
    pool = mp.Pool(processes=args.nt)
    results = pool.map(func, trj)

    pool.close()
    pool.join()

    print("Collecting data...")

    # Combine output, calculates std.
    n_length = len(results[0])
    data = [[] for _ in range(n_length)]
    for r in results:
        for n in range(n_length):
            data[n].extend(r[n])

    print("Preparing to save...")
    # n, r^2, std(r^2)
    output_data = numpy.zeros(shape=(n_length, 3))
    for n in range(n_length):
        output_data[n][0] = n
        output_data[n][1] = numpy.average(data[n])
        output_data[n][2] = numpy.std(data[n], ddof=1) / numpy.sqrt(len(data[n]))

    output = "{}msd_{}".format(args.output_prefix, "".join(args.in_file.split(".")[0:-1]))
    if args.output:
        output = args.output
    if args.csv:
        output += ".csv"
        numpy.savetxt(output, output_data)
    else:
        numpy.save(output, output_data)
    print(("Saving to {}...".format(output)))
