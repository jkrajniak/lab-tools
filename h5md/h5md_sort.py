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
import sys


def _args():
    parser = argparse.ArgumentParser('Sort H5MD file according to /id dataset')
    parser.add_argument('in_file')
    parser.add_argument('--force', action='store_true', default=False)

    return parser.parse_args()


def sort_file(h5):
    atom_groups = [ag for ag in h5['/particles'] if 'id' in h5['/particles/{}/'.format(ag)]]
    T = len(h5['/particles/{}/id/value'.format(atom_groups[0])])
    # Iterate over time frames.
    for t in xrange(T):
        sys.stdout.write('Progress: {:.2f} %\r'.format(100.0*float(t)/T))
        sys.stdout.flush()
        for ag in atom_groups:  # Iterate over atom_groups
            # Create sorted id mapping.
            ids = h5['/particles/{}/id/value'.format(ag)]
            idd = [
                x[1] for x in sorted(
                    [(p_id, col_id) for col_id, p_id in enumerate(ids[t])],
                    key=lambda y: (True, y[0]) if y[0] == -1 else (False, y[0]))
                ]
            # Iterate over group with value dataset and sort it according to particle.
            for k in h5['/particles/{}/'.format(ag)].keys():
                if 'value' in h5['/particles/{}/{}'.format(ag, k)].keys():
                    path = '/particles/{}/{}/value'.format(ag, k)
                    h5[path][t] = h5[path][t][idd]


def main():
    args = _args()
    h5 = h5py.File(args.in_file, 'r+')

    if 'sorted' in h5.attrs and not args.force:
	yn = raw_input('Warning, file was already sorted. Do you want to repeat it? (yes/no)')
        if yn == 'no':
            return False
    
    if not args.force:
        yn = raw_input('Do you want to sort file {}? (yes/no): '.format(args.in_file))
    else:
        yn = 'yes'
    if yn == 'yes':
        sort_file(h5)
        h5['/'].attrs['sorted'] = True
        h5.close()
        print('File {} sorted, closing'.format(args.in_file))


if __name__ == '__main__':
    main()
