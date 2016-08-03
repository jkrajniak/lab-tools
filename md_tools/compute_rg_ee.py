#! /usr/bin/env python
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
import mdtraj
from mdtraj.geometry import rg
import numpy as np

from md_libs import files_io

parser = argparse.ArgumentParser()
parser.add_argument('--trj', required=True)
parser.add_argument('--top', required=True)
parser.add_argument('--top_gromacs', required=True)
parser.add_argument('--out_rg', required=True)
parser.add_argument('--out_ee', required=True)
parser.add_argument('--select_atoms', required=True)
parser.add_argument('--number_of_chains', type=int, required=True)

args = parser.parse_args()

top = files_io.GROMACSTopologyFile(args.top_gromacs)
top.read()


pe10_traj = mdtraj.load(args.trj, top=args.top)
end_end_atoms_id = pe10_traj.topology.select(args.select_atoms)

# We have to calculate it for every molecule separetly 
rg_values = []
masses = np.array([top.atoms[x].mass for x in top.atoms])
for t in pe10_traj.xyz:
    rgs = rg._compute_rg_xyz(t.reshape(args.number_of_chains, t.shape[0]/args.number_of_chains, 3), masses)
    rg_values.append([np.average(rgs), np.std(rgs)])
np.savetxt(args.out_rg, rg_values)

ee_values = []
for t in pe10_traj.xyz:
    pairs = t[end_end_atoms_id.reshape(end_end_atoms_id.shape[0]/2, 2)]
    d = [(x - y) for x, y in pairs]
    d_abs = [np.sqrt(x.dot(x)) for x in d]
    ee_values.append(np.average(d_abs))

np.savetxt(args.out_ee, ee_values)
