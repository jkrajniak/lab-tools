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


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument('data_1')
    parser.add_argument('data_2')
    parser.add_argument('--min_max', help='Histogram range, format min:max')
    parser.add_argument('--bins', type=int, help='Number of bins')

    return parser.parse_args()


def test_squared_diff(hist_1, hist_2):
    """Test function: sum of squared differences."""

    diff = hist_2 - hist_1
    return np.sum(np.power(diff, 2))


def main():
    args = _args()

    data_1 = np.loadtxt(args.data_1)
    data_2 = np.loadtxt(args.data_2)

    min_bins, max_bins = map(float, args.min_max.split(':'))
    bins = np.arange(min_bins, max_bins, (max_bins-min_bins)/args.bins)

    # Create histograms with the same bins.
    histogram_1, _ = np.histogram(data_1, bins=bins)
    histogram_2, _ = np.histogram(data_2, bins=bins)

    print('f={}'.format(test_squared_diff(histogram_1, histogram_2)))


if __name__ == '__main__':
    main()
