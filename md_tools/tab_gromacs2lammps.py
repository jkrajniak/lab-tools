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
import math
import numpy as np


def _args():
    parser = argparse.ArgumentParser(
        description='Converts GROMACS tabulated potentials to LAMMPS format')
    parser.add_argument('input_file')
    parser.add_argument('output_file')
    parser.add_argument('--table_type', choices=('pair', 'bond', 'angle', 'dihedral'),
                        default='pair')
    parser.add_argument('--length_scale', default=1.0, type=float)

    return parser.parse_args()


def _pair_convert(input_f, output_f, args):
    col6 = input_f[:, 5]
    col4 = input_f[:, 3]

    r_col = input_f[1:, 0]*10.0  # Convert nm to A
    # Check if dr is the same
    dr = list(np.diff(r_col))
    if not np.isclose(dr, [dr[0]] * len(dr)).all():
        print dr.count(dr[0]), dr[0], len(dr), dr
        raise RuntimeError('The spacing in r column is not consistent')

    dr = dr[0]

    energy_col = col6
    if (col6 == 0).all():
        energy_col = col4

    energy_col = energy_col[1:]

    # Convert kJ/mol -> kcal/mol
    energy_col /= 4.184

    force_col = -np.gradient(energy_col, dr)

    print('Warning: the first row of table {} will be skipped'.format(args.input_file))

    output_f.write('{}\n'.format(args.input_file.split('.')[0]))
    output_f.write('N {} R {} {}\n\n'.format(energy_col.shape[0], np.min(r_col), np.max(r_col)))

    for i in xrange(energy_col.shape[0]):
        output_f.write('{} {} {} {}\n'.format(i+1, r_col[i], energy_col[i], force_col[i]))
    output_f.write('\n')

def _bond_convert(input_f, output_f, args):
    for line in input_f:
        l = line.strip()
        if not l.startswith('#') and not l.startswith('N') and l:
            sl = l.split()
            if len(sl) < 4:
                continue
            output_f.write('{} {} {}\n'.format(
                float(sl[1])*args.length_scale,
                float(sl[2])*4.184, float(sl[3])*41.84))


def _angle_convert(input_f, output_f, _):
    for line in input_f:
        l = line.strip()
        if not l.startswith('#') and not l.startswith('N') and l:
            sl = l.split()
            if len(sl) < 3:
                continue
            output_f.write('{} {} {}\n'.format(
                float(sl[1]), float(sl[2])*4.184,
                float(sl[3])*4.184*180.0/math.pi))


def _dihedral_convert(input_f, output_f, _):
    degrees = True
    nof = False
    data = []
    for line in input_f:
        l = line.strip()
        if not l.startswith('#') and not l.startswith('N') and l:
            sl = l.split()
            if len(sl) < 2:
                continue
            phi = float(sl[1])
            if degrees:
                phi = math.radians(phi) - math.pi
            if nof:
                data.append([phi, float(sl[2])*4.184, 0.0])
            else:
                output_f.write('{} {} {}\n'.format(
                    phi, float(sl[2])*4.184, float(sl[3])*4.184))
        elif l.startswith('N'):
            degrees = 'RADIANS' not in l
            nof = 'NOF' in l
    # Calculate force and then write a file.
    if data:
        for idx in range(0, len(data)-1):
            data[idx][2] = (data[idx+1][1] - data[idx][1])/(data[idx+1][0] - data[idx][0])
            output_f.write('{} {} {}\n'.format(*data[idx]))


def main():
    args = _args()
    input_f = np.loadtxt(args.input_file, comments=('#', '@'))
    output_f = open(args.output_file, 'w')

    print('Converting {} to {} of type {}. Warning, always convert to real lammps units'.format(
        args.input_file, args.output_file, args.table_type))
    print('The force will be recalculated.')

    table_type2func = {
        'pair': _pair_convert,
        #'bond': _bond_convert,
        #'angle': _angle_convert,
        #'dihedral': _dihedral_convert
    }

    if args.table_type not in table_type2func:
        raise RuntimeError('Table type {} currently not supported, only: {}'.format(args.table_type, table_type2func.keys()))

    table_type2func[args.table_type](input_f, output_f, args)

    output_f.close()

if __name__ == '__main__':
    main()
