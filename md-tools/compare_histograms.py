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
    parser.add_argument('data_1')
    parser.add_argument('data_2')
    parser.add_argument('--min_max', help='Histogram range, format min:max')
    parser.add_argument('--bins', type=int, help='Number of bins')
    parser.add_argument('--test_type', help='Name of test')
    parser.add_argument('--frames', type=int, default=None)

    return parser.parse_args()


def test_squared_diff(hist_1, hist_2):
    """Test function: sum of squared differences."""

    diff = hist_2 - hist_1
    return np.sum(np.power(diff, 2))


def chi_square(hist_1, hist_2):
    diff_1 = hist_1 - hist_2
    val = np.sum((np.power(diff_1, 2))/(hist_1+hist_2))
    ddof = len(hist_1)

    print('T = {}'.format(val))
    print('P(chi^2 > T) = {}'.format(chisqprob(val, ddof)))


def chi_square_shape(hist_1, hist_2):
    n1 = np.sum(hist_1)
    n2 = np.sum(hist_2)
    diff_1 = (hist_1/n1) - (hist_2/n2)
    sum_1 = (hist_1/(n1*n1)) + (hist_2/(n2*n2))
    val = np.sum((np.power(diff_1, 2))/sum_1)
    ddof = len(hist_1) - 1

    print ('T = {}'.format(val))
    print('P(chi^2 > T) = {}'.format(chisqprob(val, ddof)))


def main():
    args = _args()

    tests = {
        'squared_diff': test_squared_diff,
        'chi_square': chi_square,
        'chi_square_shape': chi_square_shape
    }
    if args.test_type not in tests:
        print('--test_type not found, available: {}'.format(tests.keys()))
        return
    print('Reading {}'.format(args.data_1))
    data_1 = np.loadtxt(args.data_1)
    print('Reading {}'.format(args.data_2))
    data_2 = np.loadtxt(args.data_2)
    print('Data read')

    if args.frames:
        data_1 = data_1[:args.frames]
        data_2 = data_2[:args.frames]

    min_bins, max_bins = map(float, args.min_max.split(':'))
    bins = np.arange(min_bins, max_bins, (max_bins-min_bins)/args.bins)

    # Create histograms with the same bins.
    histogram_1, _ = np.histogram(data_1, bins=bins, density=False)
    histogram_2, _ = np.histogram(data_2, bins=bins, density=False)

    print('Running test {}'.format(args.test_type))
    tests[args.test_type](histogram_1, histogram_2)


if __name__ == '__main__':
    main()
