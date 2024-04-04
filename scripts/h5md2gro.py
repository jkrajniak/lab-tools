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

from md_libs import files_io


def _args():
    parser = argparse.ArgumentParser("Update atom positions of input gro file from H5MD file.")
    parser.add_argument("--h5", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--input_gro", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--frame", type=int, default=-1)
    parser.add_argument("--unfolded", action="store_true", default=False)
    parser.add_argument("--valid_species")
    parser.add_argument("--scale_factor", default=1.0, type=float)
    parser.add_argument("--store_trajectory", action="store_true", default=False, help="Store whole trajectory")
    parser.add_argument("-b", help="Begin of trajectory", default=0, type=int)
    parser.add_argument("-e", help="End of trajectory", default=-1, type=int)
    parser.add_argument("--extend", action="store_true")
    return parser.parse_args()


def write_frame(args, in_gro, h5, frame, append):
    pos = h5["/particles/{}/position/value".format(args.group)][frame]
    valid_species = None
    try:
        species = h5["/particles/{}/species/value".format(args.group)][frame]
    except:
        species = h5["/particles/{}/species".format(args.group)][frame]

    if args.valid_species:
        valid_species = set(map(int, args.valid_species.split(",")))

    images = h5["/particles/{}/image/value".format(args.group)][frame]
    print(frame)
    try:
        box = np.array(h5["/particles/{}/box/edges/value".format(args.group)][frame])
    except:
        box = np.array(h5["/particles/{}/box/edges".format(args.group)])

    h5_pids = sorted([x for x in h5["/particles/{}/id/value".format(args.group)][frame] if x != -1])

    sorted(in_gro.atoms)
    ppid = 0
    max_gro_pid = max(in_gro.atoms)
    for pid, ppid in enumerate(h5_pids):
        p = pos[pid]
        s = species[pid]
        if ppid > max_gro_pid:
            in_gro.atoms[ppid] = files_io.Atom(
                atom_id=ppid,
                name="T{}".format(s),
                chain_name="XXX",
                chain_idx=ppid,
                position=p + images[pid] * box if args.unfolded else p,
            )
        if valid_species is None or species[pid] in valid_species:
            at_data = in_gro.atoms[ppid]
            if args.unfolded:
                at_data = at_data._replace(position=p + images[pid] * box)
            else:
                at_data = at_data._replace(position=p)
            if args.extend:
                at_data = at_data._replace(name="T{}".format(s))

            in_gro.atoms[ppid] = at_data
    time_frame = h5["/particles/{}/position/time".format(args.group)][frame]
    in_gro.title = "XXX molecule, t={}".format(time_frame)
    in_gro.write(args.output, force=True, append=append)


def main():
    args = _args()

    h5 = h5py.File(args.h5)
    in_gro = files_io.GROFile(args.input_gro)
    in_gro.scale_factor = args.scale_factor
    in_gro.read()

    # Get the number of frames
    pos = h5["/particles/{}/position/value".format(args.group)]
    nr_frames = pos.shape[0]
    print(("Total number of frames: {}".format(nr_frames)))
    if args.store_trajectory:
        end_frame = nr_frames
        if args.e == -1:
            end_frame = nr_frames
        elif args.e > nr_frames:
            raise RuntimeError("wrong end frame {} > {}".format(args.e, nf_frames))
        elif args.b > args.e:
            raise RuntimeError("Begin frame > end frame")
        else:
            end_frame = args.e
        frames = list(range(args.b, end_frame))
        for fr in frames:
            write_frame(args, in_gro, h5, fr, append=True)
    else:
        write_frame(args, in_gro, h5, args.frame, append=False)


if __name__ == "__main__":
    main()
