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
from md_libs import bonds
import numpy
import numpy as np
import sys


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max_tau", required=True, type=int)
    parser.add_argument("--N", required=True, type=int, help="Number of molecules")
    parser.add_argument("--interactive", action="store_true", default=False)
    parser.add_argument("--begin", default=0, type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--output")
    parser.add_argument("--csv", action="store_true")
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--every-frame", dest="step", default=1, type=int)
    parser.add_argument("--prefix", default="", type=str)
    parser.add_argument("--normalize_c0", action="store_true", default=False)
    parser.add_argument("in_file")

    return parser.parse_args()


def calculate_end_end_acf(trajectory, max_tau):
    """Calculating autocorrelation function of end-to-end vector
    Args:
        trajectory: numpy array with trajectory.
        max_tau: Maximum tau.
    Returns:
        Returns ACF and ACF error.
    """
    N = len(trajectory[0])  # number of molecules
    max_t = len(trajectory)

    acf = np.zeros(max_tau)
    acf_error = np.zeros(max_tau)
    sumdist = 0.0
    print("Frames: 0..", max_t)
    print("Max tau:", max_tau)
    # Go over time.
    for m in range(max_tau):
        sys.stdout.write("Tau: {}\r".format(m))
        sys.stdout.flush()
        intermediate_results = []
        for n in range(0, max_t - m):
            sumdist = 0.0
            for i in range(N):
                a1 = trajectory[n][i]
                a2 = trajectory[n + m][i]
                sumdist += np.dot(a1 / np.linalg.norm(a1), a2 / np.linalg.norm(a2))
            intermediate_results.append(sumdist / N)
        acf[m] = np.average(intermediate_results)
        acf_error[m] = np.std(intermediate_results, ddof=1)  # /np.sqrt(max_t-m)
    return acf, acf_error


def main():
    args = _args()
    data = numpy.load(args.in_file)

    max_N = min([len(x) for x in data])
    if max_N < args.N:
        print(("max_N: {}, args.N: {}".format(max_N, args.N)))
    else:
        pass

    if not args.end:
        args.end = len(data)
        print(args.begin, args.end)

    raw_data = data[args.begin : args.end : args.step]

    print("Normalize vectors")
    data = np.zeros(raw_data.shape)
    norms = np.apply_along_axis(np.linalg.norm, 2, raw_data)
    print(data.shape, norms.shape)
    for i in range(data.shape[0]):
        data[i] = (raw_data[i].T / norms[i]).T

    args.dt *= args.step
    max_tau = args.max_tau / args.step
    print(("Traj length: {}".format(len(data))))
    print(("dt: {}".format(args.dt)))
    print(("max_tau: {}".format(max_tau)))

    print("Calculating ACF...")
    # acf1, acf_errors = calculate_end_end_acf(data, max_tau)
    acf1, acf_errors = bonds.calculate_end_end_acf(data, max_tau)

    if args.normalize_c0:
        print("Normalize data by C(0)")
        output_data = acf1 / acf1[0]
    else:
        output_data = acf1
    time_column = numpy.arange(len(output_data))  # *args.dt
    output = "{}acf_{}".format(args.prefix, args.in_file)
    if args.output:
        output = args.output
    print(("Saving to {}...".format(output)))
    if args.csv:
        output += ".csv"
        numpy.savetxt(output, numpy.column_stack((time_column, output_data, acf_errors)))
    else:
        numpy.save(output, numpy.column_stack((time_column, output_data, acf_errors)))

    if args.interactive:
        import IPython

        IPython.embed()


if __name__ == "__main__":
    main()
