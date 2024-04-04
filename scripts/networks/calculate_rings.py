#!/usr/bin/env python3
"""
Copyright (C) 2018 Jakub Krajniak <jkrajniak@gmail.com>

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
import h5py
import networkx as nx
import numpy as np
import multiprocessing as mp


DYNAMIC_BONDS = "/connectivity/chem_bonds_0"
STATIC_BONDS = "/connectivity/bonds_0"


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5")
    parser.add_argument("output")
    parser.add_argument("--nt", default=None, type=int)

    return parser.parse_args()


def compute_frame(h5filename, time_idx):
    with h5py.File(h5filename, "r", driver="stdio") as h5:
        pid2res_id = {i: v for i, v in enumerate(h5["/particles/atoms/res_id/value"][time_idx], 1) if v != -1}
        static_bonds = np.array(h5[STATIC_BONDS])
        static_bonds = static_bonds[np.where(static_bonds != [-1, -1])[0]]
        dynamic_bonds = np.array(h5[DYNAMIC_BONDS + "/value"][time_idx])
        dynamic_bonds = dynamic_bonds[np.where(dynamic_bonds != [-1, -1])[0]]
        graph_bonds = {
            tuple(sorted(map(pid2res_id.get, x)))  ## Build super CG model
            for x in dynamic_bonds
        }
        # for x in np.concatenate((static_bonds, dynamic_bonds))}
        graph_bonds = {b for b in graph_bonds if b[0] != b[1]}
        graph = nx.Graph(
            incoming_graph_data=list(graph_bonds),
            **{"conversion": dynamic_bonds.shape[0] / 2000.0, "num_bonds": dynamic_bonds.shape[0]},
        )
        print(f"Frame idx {time_idx}")
        return compute_simple_loops(graph)


def compute_simple_loops(g, g_sub=None):
    """Gets the list of loops without branches."""
    if g_sub is None:
        g_sub = list(nx.connected_component_subgraphs(g))

    loop_size = []
    for sub_g in g_sub:
        if sub_g.number_of_nodes() == 1:
            continue
        cycle_set = [set(sub_g.subgraph(k).edges()) for k in nx.cycle_basis(sub_g)]
        if len(cycle_set) > 1:
            for i in range(len(cycle_set)):
                for j in range(i + 1, len(cycle_set)):
                    cy_01 = list(set.symmetric_difference(cycle_set[i], cycle_set[j]))
                    g_01 = nx.Graph()
                    g_01.add_edges_from(cy_01)
                    for g_01s in nx.connected_component_subgraphs(g_01):
                        loop_size.append(g_01s.number_of_nodes())
        elif len(cycle_set) == 1:
            loop_size.append(len({n for e in cycle_set[0] for n in e}))

    if loop_size:
        return (loop_size, np.average(loop_size), len(loop_size), g.graph["conversion"], len(g_sub), g)
    return ([], 0.0, 0, g.graph["conversion"], len(g_sub), g)


def main():
    args = _args()
    with h5py.File(args.h5, "r") as h5:
        cl = h5[DYNAMIC_BONDS + "/value"]
        frames = range(cl.shape[0])

    compute_frame_ = functools.partial(compute_frame, args.h5)

    p = mp.Pool(args.nt)
    loops = p.map(compute_frame_, frames)
    avg_loop_size = [[x[3], x[1], x[2], x[4]] for x in loops]
    np.savetxt(args.output, avg_loop_size, header="p avg_loop_size num_loops num_components")


if __name__ == "__main__":
    main()
