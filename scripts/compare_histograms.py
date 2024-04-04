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
import numpy as np
from scipy.stats import chisqprob


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_1")
    parser.add_argument("data_2")
    parser.add_argument("--raw_data", action="store_true", default=False, help="data_1 and data_2 are a raw data")
    parser.add_argument("--min_max", help="Histogram range, format min:max")
    parser.add_argument("--bins", type=int, help="Number of bins")
    parser.add_argument("--test_type", help="Name of test")
    parser.add_argument("--frames", type=int, default=None)

    return parser.parse_args()


def test_squared_diff(hist_1, hist_2):
    """Test function: sum of squared differences."""

    diff = hist_2 - hist_1
    print(("T = {}".format(np.nansum(np.power(diff, 2)))))


def chi_square(hist_1, hist_2):
    diff_1 = hist_1 - hist_2
    val = np.nansum((np.power(diff_1, 2)) / (hist_1 + hist_2))
    ddof = len(hist_1)

    print(("T = {}".format(val)))
    print(("P(chi^2 > T) = {}".format(chisqprob(val, ddof))))


def chi_square_shape(hist_1, hist_2):
    n1 = np.sum(hist_1)
    n2 = np.sum(hist_2)
    diff_1 = (hist_1 / n1) - (hist_2 / n2)
    sum_1 = (hist_1 / (n1 * n1)) + (hist_2 / (n2 * n2))
    val = np.nansum((np.power(diff_1, 2)) / sum_1)
    ddof = len(hist_1) - 1

    print(("T = {}".format(val)))
    print(("P(chi^2 > T) = {}".format(chisqprob(val, ddof))))


def ecf(hist_1, bins):
    """Computes empirical cumulative distribution."""
    n = sum(hist_1[:, 1])
    return np.array([np.sum([y for x, y in hist_1 if x < xl]) for xl in bins]) / n


def kolmogorov_smirnov(hist_1, hist_2):
    """Kolmogorov-Smirnov"""
    n1 = len(hist_1)
    n2 = len(hist_2)
    ecf_1 = ecf(hist_1, hist_1[:, 0])
    ecf_2 = ecf(hist_2, hist_1[:, 0])
    distance = np.max(ecf_1 - ecf_2)
    pref = np.sqrt((n1 * n2) / (n1 + n2))
    c_a = 1.36
    accept_null = distance < c_a / np.sqrt(n1)
    print(("T = {}".format(distance)))
    print(("pref = {}".format(pref)))
    print(("Accept null = {}".format(accept_null)))


def main():
    args = _args()

    tests = {
        "squared_diff": test_squared_diff,
        "chi_square": chi_square,
        "chi_square_shape": chi_square_shape,
        "ks": kolmogorov_smirnov,
    }
    if args.test_type not in tests:
        print(("--test_type not found, available: {}".format(list(tests.keys()))))
        return

    # Create histograms with the same bins.
    if args.raw_data:
        print(("Reading {}".format(args.data_1)))
        data_1 = np.loadtxt(args.data_1)
        print(("Reading {}".format(args.data_2)))
        data_2 = np.loadtxt(args.data_2)
        print("Data read")

        if args.frames and args.raw_data:
            data_1 = data_1[: args.frames]
            data_2 = data_2[: args.frames]

        min_bins, max_bins = list(map(float, args.min_max.split(":")))
        bins = np.arange(min_bins, max_bins, (max_bins - min_bins) / args.bins)
        histogram_1, _ = np.histogram(data_1, bins=bins, density=False)
        histogram_2, _ = np.histogram(data_2, bins=bins, density=False)
    else:
        histogram_1 = np.loadtxt(args.data_1, usecols=(0, 1))
        histogram_2 = np.loadtxt(args.data_2, usecols=(0, 1))

    print(("Running test {}".format(args.test_type)))
    tests[args.test_type](histogram_1, histogram_2)


if __name__ == "__main__":
    main()
