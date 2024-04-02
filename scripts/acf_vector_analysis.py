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
import datetime
import numpy as np
try:
    import MDAnalysis
    SUPPORT_GROMACS=True
except ImportError:
    SUPPORT_GROMACS=False
import time

import h5py

from md_libs import bonds as bond_libs
from md_libs import files_io
from mpi4py import MPI

size = MPI.COMM_WORLD.size
rank = MPI.COMM_WORLD.rank


class ListAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        setattr(
            namespace,
            option_string.replace('-', ''),
            [list(map(int, x.split('-'))) for x in values.split(',')]
        )


def _args():
    parser = argparse.ArgumentParser('Analyze bond distance')
    parser.add_argument('--trj', help='Trajectory file', required=True)
    parser.add_argument('--top', help='Topology file')
    parser.add_argument('--molecules', help='Number of molecules', required=True, type=int)
    parser.add_argument('--N', help='Number of atoms in molecule', required=True, type=int)
    parser.add_argument('--group', help='Atoms group (only for H5MD)', required=False)
    parser.add_argument('--vector', help='Vector, in format 1-12', action=ListAction, default=[],
                        required=True)
    parser.add_argument('--box_left', help='Bottom front left point of the box, format (x, y, z)')
    parser.add_argument('--box_right', help='Top rear right point of the box, format (x, y, z)')
    parser.add_argument('--begin', help='Begin frame', default=0, type=int)
    parser.add_argument('--end', help='End frame', default=-1, type=int)
    parser.add_argument('--prefix', help='Prefix')
    args = parser.parse_args()

    return args


def replicate_list(input_list, mol, N, shift=-1, cmplx=False):
    rr = [[i+x*N+shift for i in v] for x in range(mol) for v in input_list]
    if cmplx:
        return rr
    else:
        return [x for v in rr for x in v]


def _gromacs_processing(args, q_vector, start_stop_list, pos_cons):
    pid = MPI.COMM_WORLD.rank
    if pid >= len(start_stop_list):
        return [], [], []
    start_stop = start_stop_list[pid]
    vectors = []

    process_tuples = [
        [vectors, q_vector, 2, bond_libs.calculate_bond_vec]
    ]

    trj = MDAnalysis.coordinates.reader(args.trj)
    box = np.array(trj.ts.dimensions[:3])
    half_box = 0.5*box

    if args.end == -1:
        args.end = trj.n_frames

    MDAnalysis.core.flags['use_periodic_selections'] = True
    MDAnalysis.core.flags['use_KDTree_routines'] = False

    print('%d processing %s, pos_cons %s' % (pid, start_stop, pos_cons))
    frame_idx = start_stop[0]
    for t in trj[start_stop[0]:start_stop[1]]:
        print('Frame %d' % frame_idx)
        for output, atom_ids, tuple_size, functor in process_tuples:
            if atom_ids is None or atom_ids.size == 0:
                continue
            atoms = np.array(t)[atom_ids]
            # filter with pos_cons
            atom_tuples = []
            for d_tuple in atoms:
                valid = False
                for a in d_tuple:
                    if ((a[0] >= pos_cons[0][0] and a[0] <= pos_cons[1][0]) and
                        (a[1] >= pos_cons[0][1] and a[1] <= pos_cons[1][1]) and
                        (a[2] >= pos_cons[0][2] and a[2] <= pos_cons[1][2])):
                        valid = True
                    else:
                        valid = False
                        break
                if valid:
                    atom_tuples.append([np.array(x) for x in d_tuple])
                else:
                    atom_tuples.append([None for x in d_tuple])
            output.append(functor(np.array(atom_tuples), box, half_box))
        frame_idx += 1
    return vectors


def _h5_processing(args, q_vector, start_stop_list, pos_cons):
    nt = MPI.COMM_WORLD.size
    at_group = args.group
    if nt > 1:
        pid = MPI.COMM_WORLD.rank
        if pid >= len(start_stop_list):
            return []
        start_stop = start_stop_list[pid]
        h5file = h5py.File(args.trj, 'r', driver='mpio', comm=MPI.COMM_WORLD)
    else:
        h5file = h5py.File(args.trj, 'r')
        start_stop = start_stop_list[0]

    _, box, trj, _ = files_io.prepare_h5md(h5file, at_group, start_stop[0], start_stop[1])    

    half_box = 0.5*box

    vectors = []
    

    process_tuples = [
        [vectors, q_vector, bond_libs.calculate_bond_vec],
    ]

    frame_idx = start_stop[0]
    for frame in trj:
        print(('Frame {}'.format(frame_idx)))
        for output, atom_ids, functor in process_tuples:
            if atom_ids is None or atom_ids.size == 0:
                continue
            atoms = frame[atom_ids]
            atom_tuples = []
            for d_tuple in atoms:
                atom_tuples.append([np.array(x) for x in d_tuple])
            output.append(functor(np.array(atom_tuples), box, half_box))
        frame_idx += 1
    return vectors


def _get_info(args):
    """Return information about box and trajectory."""

    if (args.trj.endswith('trr') or args.trj.endswith('xtc')) and SUPPORT_GROMACS:
        trj = MDAnalysis.coordinates.reader(args.trj)
        box = np.array(trj.ts.dimensions[:3])
        return box, trj.n_frames, True
    elif args.trj.endswith('h5'):
        h5file = h5py.File(args.trj, 'r')
        trj = h5file['particles/{}/position/value'.format(args.group)]
        box = np.array(h5file['particles/{}/box/edges'.format(args.group)])
        return box, len(trj), True
    else:
        raise RuntimeError('Wrong trajectory')



def save_data(name, array, definitions, filename_prefix, file_template):
    for idx, definition in enumerate(definitions):
        filename = '%s%s_%s.dat' % (
            filename_prefix, name, '_'.join(map(str, definition)))
        print(('Saving file {}.npy'.format(filename)))
        np.save(filename, np.array([
            np.array([x for x in arr[idx::len(definitions)] if x is not None])
            for arr in array
        ]))


def main_local():
    args = _args()
    time0 = time.time()

    box, numframes, cmplx = _get_info(args)
    if args.end == -1:
        args.end = numframes

    print(('Box: {}'.format(box)))
    print(('Frames: {}'.format(numframes)))

    if args.box_left and args.box_right:
        x0, y0, z0 = list(map(float, args.box_left.split(',')))
        x1, y1, z1 = list(map(float, args.box_right.split(',')))
        pos_cons = ((x0, y0, z0), (x1, y1, z1))
    else:
        pos_cons = ((0.0, 0.0, 0.0), tuple(box))

    # Prepare for MPI run
    nt = MPI.COMM_WORLD.size
    if nt > 1:
        num_frames = args.end - args.begin
        per_process = int(round(num_frames / float(nt)))
        print('per_process', per_process)
        frame_ranges = list(range(args.begin, args.end, per_process))
        frame_range_list = [
            (frame_ranges[i], frame_ranges[i+1]) for i in range(len(frame_ranges)-1)
        ]
        frame_range_list.append((frame_ranges[-1], args.end+1))  # Because it is right open range
    else:
        frame_range_list = [(args.begin, args.end)]

    print(('Frame range: {}'.format(frame_range_list)))
    # Q
    q_vector = np.array(replicate_list(args.vector, args.molecules, args.N, cmplx=cmplx))

    # PMI run
    import pmi
    pmi.setup()
    pmi.execfile_(__file__)
    if (args.trj.endswith('trr') or args.trj.endswith('xtc')) and SUPPORT_GROMACS:
        data = pmi.invoke(
            '_gromacs_processing',
            args, q_vector, frame_range_list, pos_cons)
    elif args.trj.endswith('h5'):
        data = pmi.invoke(
            '_h5_processing',
            args, q_vector, frame_range_list, pos_cons)
    else:
        raise RuntimeError('Wrong trajectory file')

    # Collect datas
    vectors = []
    for node_data in data:
        vectors.extend([v for v in node_data])

    filename_prefix = '' if not args.prefix else args.prefix + '_'
    file_template = '# Date: %s\n# Filename: %s\n' % (datetime.datetime.today(), args.trj)
    print('Saving data...')
    save_data('vector', vectors, args.vector, filename_prefix, file_template)
    print(('Processing time: {}s with {} CPUs'.format(time.time() - time0, nt)))

if __name__ != 'pmi':
    main_local()
