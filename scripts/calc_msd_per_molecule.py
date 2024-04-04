#!/usr/bin/env python
"""
Copyright (C) 2017 Jakub Krajniak <jkrajniak@gmail.com>

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
import functools
import numpy as np
import os
import subprocess
import multiprocessing as mp

from md_libs import files_io


def calc_msd(traj_file, output_prefix, msd_command, args):
    atom_ids, mol_id = args
    current = mp.current_process()
    rank = current._identity[0]
    index_file = "tmp_water_{}.ndx".format(rank)
    with open(index_file, "w") as iw:
        iw.write("[ H2O ]\n")
        iw.write(" ".join(map(str, atom_ids)))
        iw.write("\n")

    DEVNULL = open(os.devnull, "wb")
    cmds = msd_command.split() + ["-f", traj_file, "-n", index_file, "-o", "{}msd_{}.xvg".format(output_prefix, mol_id)]
    msd_cmd = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=DEVNULL)
    data = msd_cmd.stdout.readlines()
    ret = [mol_id, -1, -1]
    if data:
        msd_part = data[-1].split()
        msd_data = list(map(float, [msd_part[2], msd_part[4].replace(")", "")]))
        ret = [mol_id, msd_data[0], msd_data[1]]
    print(ret)
    os.remove(index_file)
    return ret


def _args():
    parser = argparse.ArgumentParser("Calculate MSD for individual molecules")
    parser.add_argument("--in_top", help=".top file", required=True)
    parser.add_argument("--trj", help="Trajectory file", required=True)
    parser.add_argument("--nt", default=4, type=int, help="Number of processes")
    parser.add_argument("--mol_name", help="Water molecule name")
    parser.add_argument("--mol_range", help="Molecule range, e.g. start:stop")
    parser.add_argument("--mol_size", default=3, type=int, help="Size of water molecule")
    parser.add_argument("--msd_command", default="gmx_mpi msd")
    parser.add_argument("--output", default="msd_data", help="Output file")
    parser.add_argument("--output_prefix", default="", help="Prefix for files")

    return parser


def main():
    args = _args().parse_args()

    top = files_io.GROMACSTopologyFile(args.in_top)
    top.read()

    water_molecules = [x for x in sorted(top.atoms) if top.atoms[x].chain_name == args.mol_name]

    if args.mol_name:
        atom_ids = [
            (water_molecules[i : i + args.mol_size], mol_id)
            for mol_id, i in enumerate(range(0, len(water_molecules), args.mol_size))
        ]
    elif args.mol_range:
        start, stop = list(map(int, args.mol_range.split(":")))
        mol_ids = list(range(start, stop))
        water_molecules = [x for x in sorted(top.atoms) if top.atoms[x].chain_idx in mol_ids]
        atom_ids = [
            (water_molecules[i : i + args.mol_size], mol_id)
            for mol_id, i in enumerate(range(0, len(water_molecules), args.mol_size))
        ]
    else:
        raise RuntimeError("Please specify either mol_name or mol_range.")

    print(("Selected atoms {}".format(len(atom_ids))))

    p = mp.Pool(args.nt)

    _calc_msd = functools.partial(calc_msd, args.trj, args.output_prefix, args.msd_command)
    mol_msd_data = p.map(_calc_msd, atom_ids)

    save_file = "{}{}".format(args.output_prefix, args.output)
    np.savetxt(save_file, mol_msd_data, header="mol_id D +/-")
    print(("Saved {}".format(save_file)))


if __name__ == "__main__":
    main()
