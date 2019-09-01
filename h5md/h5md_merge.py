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
import numpy as np
import shutil

parser = argparse.ArgumentParser(
    'Stack input file on top of output file (only groups in /particles)')
parser.add_argument('input_file', help='Input H5MD file')
parser.add_argument('output_file', help='Output (appendable) H5MD file')
parser.add_argument('group_name', help='Particles group')
parser.add_argument('--backup_output', action='store_false', default=True)
parser.add_argument('--renumber_steps', action='store_true',
                    help='Renumber /step dataset (with the same spacing!')
parser.add_argument('--renumber_time', action='store_true',
                    help='Renumber /time dataset (with the same spacing!')

args = parser.parse_args()

# Merge only particles
if args.backup_output:
    print('Making backup {file}->{file}.bak'.format(file=args.output_file))
    shutil.copy(args.output_file, '{}.bak'.format(args.output_file))

input_file = h5py.File(args.input_file, 'r')
output_file = h5py.File(args.output_file, 'r+')
group_name = args.group_name

old_particles = output_file['/particles/{}/'.format(group_name)]
new_particles = input_file['/particles/{}/'.format(group_name)]


def merge_datasets(path):
    input_dataset = new_particles[path]
    output_dataset = old_particles[path]
    old_shape = output_dataset.shape
    output_dataset.resize(output_dataset.shape[0] + input_dataset.shape[0], 0)
    output_dataset[old_shape[0]:] = input_dataset
    return output_dataset


def renumber_dataset(ds):
    if len(ds.shape) > 1:
        return
    spacing = ds[3] - ds[2]
    ds_length = ds.shape[0]
    start_element = ds[0]
    ds[:] = np.arange(start_element, start_element+ds_length*spacing, spacing, dtype=ds.dtype)
    assert np.all(np.diff(ds) == spacing)


def merge_group(name, object):
    if isinstance(object, h5py._hl.group.Group):
        return
    print('Merging dataset {}'.format(object.name))
    result_dataset = merge_datasets(object.name)
    if args.renumber_steps and name.endswith('step'):
        print('Renumbering {}'.format(object.name))
        renumber_dataset(result_dataset)
    if args.renumber_time and name.endswith('time'):
        print('Renumbering {}'.format(object.name))
        renumber_dataset(result_dataset)

old_particles.visititems(merge_group)
