#!/usr/bin/env python
"""
Copyright (C) 2014 Jakub Krajniak <jkrajniak@gmail.com>

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
from collections import namedtuple
import datetime
import numpy as np
try:
    import MDAnalysis
    SUPPORT_GROMACS=True
except:
    SUPPORT_GROMACS=False
import time

import h5py

from libs import bonds as bond_libs
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
    parser.add_argument('--bonds', help='List of bonds, ex. 1-2,2-3,3-4,1-4', action=ListAction,
                        default=[])
    parser.add_argument('--angles', help='List of angles, ex. 1-2-3,2-3-4', action=ListAction,
                        default=[])
    parser.add_argument('--torsions', help='List of torsions, ex. 1-2-3-4', action=ListAction,
                        default=[])
    parser.add_argument('--box_left', help='Bottom front left point of the box, format (x, y, z)')
    parser.add_argument('--box_right', help='Top rear right point of the box, format (x, y, z)')
    parser.add_argument('--begin', help='Begin frame', default=0, type=int)
    parser.add_argument('--end', help='End frame', default=-1, type=int)
    parser.add_argument('--prefix', help='Prefix')
    parser.add_argument('--timeseries', action='store_true', default=False)
    parser.add_argument('--scalling', default=1.0, type=float)
    args = parser.parse_args()

    return args


def replicate_list(input_list, mol, N, shift=-1, cmplx=False):
    rr = [[i+x*N+shift for i in v] for x in range(mol) for v in input_list]
    if cmplx:
        return rr
    else:
        return [x for v in rr for x in v]


def _gromacs_processing(args, q_bonds, q_angles, q_torsions, start_stop_list, pos_cons, box):
    pid = MPI.COMM_WORLD.rank
    if pid >= len(start_stop_list):
        return [], [], []
    start_stop = start_stop_list[pid]
    bonds = []
    angles = []
    torsions = []

    process_tuples = [
        [bonds, q_bonds, 2, bond_libs.calculate_bond],
        [angles, q_angles, 3, bond_libs.calculate_angle],
        [torsions, q_torsions, 4, bond_libs.calculate_dihedral]
    ]
    
    trajectory = MDAnalysis.coordinates.reader(args.trj)
    half_box = 0.5*box

    if args.end == -1:
        args.end = trajectory.numframes

    MDAnalysis.core.flags['use_periodic_selections'] = True
    MDAnalysis.core.flags['use_KDTree_routines'] = False

    print('%d processing %s, pos_cons %s' % (pid, start_stop, pos_cons))
    frame_idx = start_stop[0]
    for t in trajectory[start_stop[0]:start_stop[1]]:
        print('Frame %d' % frame_idx)
        for output, atom_ids, tuple_size, functor in process_tuples:
            if atom_ids is None or atom_ids.size == 0:
                continue
            atoms = np.array(t)[atom_ids]
            atom_tuples = []
            for d_tuple in atoms:
                valid = False
                for a in d_tuple:
                    a *= args.scalling
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

            atom_tuples = np.array(atom_tuples)
            if args.timeseries:
                output.append(functor(atom_tuples, box, half_box))
            else:
                output.extend(functor(atom_tuples, box, half_box))

        frame_idx += 1

    return bonds, angles, torsions


def _h5_processing(args, q_bonds, q_angles, q_torsions, start_stop_list, pos_cons, box):
    nt = MPI.COMM_WORLD.size
    at_group = args.group
    if nt > 1:
        pid = MPI.COMM_WORLD.rank
        if pid >= len(start_stop_list):
            return [], [], []
        start_stop = start_stop_list[pid]
        h5file = h5py.File(args.trj, 'r', driver='mpio', comm=MPI.COMM_WORLD)
    else:
        h5file = h5py.File(args.trj, 'r')
        start_stop = start_stop_list[0]

    bonds, angles, torsions = [], [], []

    trj = h5file['particles/{}/position/value'.format(at_group)][start_stop[0]:start_stop[1]]
    half_box = 0.5*box

    bonds, angles, torsions = [], [], []

    process_tuples = [
        [bonds, q_bonds, bond_libs.calculate_bond],
        [angles, q_angles, bond_libs.calculate_angle],
        [torsions, q_torsions, bond_libs.calculate_dihedral]
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
                valid = False
                for a in d_tuple:
                    a *= args.scalling
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
            atom_tuples = np.array(atom_tuples)
            if args.timeseries:
                output.append(functor(atom_tuples, box, half_box))
            else:
                output.extend(functor(atom_tuples, box, half_box))
        frame_idx += 1

    return bonds, angles, torsions


def _get_info(args):
    """Return information about box and trajectory."""
    return_tuple = namedtuple('Return', ['box', 'num_frames', 'cmplx', 'is_gromacs'])

    if (args.trj.endswith('trr') or args.trj.endswith('xtc') or args.trj.endswith('gro')) and SUPPORT_GROMACS:
        trj = MDAnalysis.coordinates.reader(args.trj)
        box = np.array(trj.ts.dimensions[:3])*args.scalling
        return return_tuple(box, trj.numframes, True, True)
    elif args.trj.endswith('h5'):
        h5file = h5py.File(args.trj, 'r')
        trj = h5file['particles/{}/position/value'.format(args.group)]
        box = np.array(h5file['particles/{}/box/edges/value'.format(args.group)][-1])*args.scalling
        return return_tuple(box, len(trj), True, False)
    else:
        raise RuntimeError('Wrong trajectory')


def save_data(name, array, definitions, filename_prefix, file_template, timeserie):
    for idx, definition in enumerate(definitions):
        filename = '%s%s_%s.dat' % (
            filename_prefix, name, '_'.join(map(str, definition)))
        print(('Saving file {}'.format(filename)))
        if timeserie:
            np.save(filename, np.array([
                np.array([x for x in arr[idx::len(definitions)] if x is not None])
                for arr in array
            ]))
        else:
            with open(filename, 'w') as fout:
                fout.write(file_template)
                fout.write('# %s %s\n' % (name, definition))
                fout.write('\n')
                np.savetxt(fout, [x for x in array[idx::len(definitions)] if x is not None])


def main_local():
    args = _args()
    time0 = time.time()

    box, numframes, cmplx, is_gromacs = _get_info(args)
    if args.end == -1:
        args.end = numframes

    print(('Box: {}'.format(box)))
    print(('Frames: {}'.format(numframes)))
    print(('Scalling factor: {}'.format(args.scalling)))

    if args.box_left and args.box_right:
        x0, y0, z0 = [float(x)*args.scalling for x in args.box_left.split(',')]
        x1, y1, z1 = [float(x)*args.scalling for x in args.box_right.split(',')]
        pos_cons = ((x0, y0, z0), (x1, y1, z1))
    else:
        pos_cons = ((0.0, 0.0, 0.0), tuple(box))
    print(('Position constraints: {}'.format(pos_cons)))

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
    q_bonds = np.array(replicate_list(args.bonds, args.molecules, args.N, cmplx=cmplx))
    q_angles = np.array(replicate_list(args.angles, args.molecules, args.N, cmplx=cmplx))
    q_torsions = np.array(replicate_list(args.torsions, args.molecules, args.N, cmplx=cmplx))

    # PMI run
    import pmi
    pmi.setup()
    pmi.execfile_(__file__)
    if is_gromacs:
        data = pmi.invoke(
            '_gromacs_processing',
            args, q_bonds, q_angles, q_torsions, frame_range_list, pos_cons, box)
    else:
        data = pmi.invoke(
            '_h5_processing',
            args, q_bonds, q_angles, q_torsions, frame_range_list, pos_cons, box)

    # Collect datas
    bonds = []
    angles = []
    torsions = []
    for node_data in data:
        bonds.extend([v for v in node_data[0]])
        angles.extend([v for v in node_data[1]])
        torsions.extend([v for v in node_data[2]])

    filename_prefix = '' if not args.prefix else args.prefix + '_'
    file_template = '# Date: %s\n# Filename: %s\n' % (datetime.datetime.today(), args.trj)
    print('Saving data...')
    save_data('bond', bonds, args.bonds, filename_prefix, file_template, args.timeseries)
    save_data('angle', angles, args.angles, filename_prefix, file_template, args.timeseries)
    save_data('torsion', torsions, args.torsions, filename_prefix, file_template, args.timeseries)
    print(('Processing time: {}s with {} CPUs'.format(time.time() - time0, nt)))

if __name__ != 'pmi':
    main_local()
