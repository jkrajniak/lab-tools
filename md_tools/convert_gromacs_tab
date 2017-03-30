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
    parser = argparse.ArgumentParser('Convert GROMACS tabulated potential')
    parser.add_argument('in_tab', help='Input table')
    parser.add_argument('out_tab', help='Output table')
    gr = parser.add_mutually_exclusive_group(required=True)
    gr.add_argument('--c12c6', help='Move data from 6-th column to 4-th column',
                    action='store_true')
    gr.add_argument('--c6c12', help='Move data from 4-th columnt to 6-th column',
                    action='store_true')
    parser.add_argument('--skip_lines', default=0, help='Number of lines to skip in input file',
                        type=int)
    parser.add_argument('--comment_chars', default='#',
                        help='Skip lines that starts with characters')

    return parser.parse_args()


def main():
    args = _args()

    in_tab = np.loadtxt(args.in_tab, comments=list(args.comment_chars), skiprows=args.skip_lines)

    if args.c12c6:
        in_tab[:, (3, 5)] = in_tab[:, (5, 3)]
        in_tab[:, (4, 6)] = in_tab[:, (6, 4)]
    elif args.c6c12:
        in_tab[:, (5, 3)] = in_tab[:, (3, 5)]
        in_tab[:, (6, 4)] = in_tab[:, (4, 6)]

    np.savetxt(args.out_tab, in_tab, fmt='%.10e', header='Converted from {}'.format(args.in_tab))
    print('File saved {}'.format(args.out_tab))

if __name__ == '__main__':
    main()

