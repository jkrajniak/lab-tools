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
import functools
import re
from md_libs import files_io
import networkx as nx

import multiprocessing


def _args():
    parser = argparse.ArgumentParser("Generate exclusion lists")
    parser.add_argument("in_top", help="Input topology to process")
    parser.add_argument("out_list", help="Output exclusion list in the form of pairs")
    parser.add_argument("--nrexcl", type=int, default=3)
    parser.add_argument(
        "--select_atoms",
        default=None,
        help=(
            "Generate exclusion list for subgroup of atoms. This is a comma separated"
            " list of regular expressions (compiled by Python re module)"
        ),
    )

    parser.add_argument("--nt", default=None, help="Number of cores", type=int)
    parser.add_argument("--append", default=False, action="store_true", help="Append to the file")

    return parser.parse_args()


def get_all_paths(g_input, nrexcl, i_node):
    node_list = sorted(g_input.node)
    paths = []
    for j in range(g_input.number_of_nodes()):
        if i_node != j:
            p = nx.all_simple_paths(g_input, node_list[i_node], node_list[j], nrexcl)
            paths.extend([tuple(sorted((x[0], x[-1]))) for x in p if x[0] != x[-1]])
    return paths


def replicate_list(input_list, mol, N, shift=-1, cmplx=False):
    rr = [[i + x * N + shift for i in v] for x in range(mol) for v in input_list]
    if cmplx:
        return rr
    else:
        return [x for v in rr for x in v]


def main():
    args = _args()
    input_top = files_io.GROMACSTopologyFile(args.in_top)
    input_top.read()

    print(("Generate exclusion list with nrexcl={}".format(args.nrexcl)))

    selected_atom_ids = list(input_top.atoms.keys())
    if args.select_atoms:
        selected_atom_ids = []
        select_list = list(map(re.compile, args.select_atoms.split(",")))
        for at_id, at_data in list(input_top.atoms.items()):
            val = [re.match(m, at_data.name) for m in select_list]
            if any(val):
                selected_atom_ids.append(at_id)
    print(("Selected atom symbols: {}".format(",".join({input_top.atoms[x].name for x in selected_atom_ids}))))

    # Build a graph
    g = nx.Graph()
    for b1, b2 in list(input_top.bonds.keys()):
        if b1 in selected_atom_ids and b2 in selected_atom_ids:
            g.add_edge(b1, b2)

    for b1, b2 in list(input_top.cross_bonds.keys()):
        if b1 in selected_atom_ids and b2 in selected_atom_ids:
            g.add_edge(b1, b2)

    print(("Generated graph, number of edges {}, nodes {}".format(g.number_of_edges(), g.number_of_nodes())))
    all_paths = []
    connected_subgraphs = nx.connected_component_subgraphs(g)
    pool = multiprocessing.Pool(args.nt)
    for sub_g in connected_subgraphs:
        sorted(sub_g.node)
        path_gen = functools.partial(get_all_paths, sub_g, args.nrexcl)
        ans = pool.map(path_gen, list(range(sub_g.number_of_nodes())))
        all_paths.extend([x for l in ans for x in l])

    all_paths = sorted(set(all_paths))
    print(("Generate {} exclusions".format(len(all_paths))))

    # If topology contains description of nmols then we have to replicate
    # the exclusion lists.
    replicated_paths = []
    for molecule in input_top.molecules:
        at_id = [x for x, d in list(input_top.atoms.items()) if d.chain_name == molecule["name"]]
        nmols = int(molecule["mol"])
        n_at = len(at_id)
        print(("Process molecule: {}, atoms {}".format(molecule, len(at_id))))
        replicated_paths.extend(
            replicate_list([x for x in all_paths if x[0] in at_id and x[1] in at_id], nmols, n_at, cmplx=True, shift=0)
        )

    print(("Generate {} exclusions".format(len(replicated_paths))))
    with open(args.out_list, "a" if args.append else "w") as output_file:
        for p in replicated_paths:
            output_file.write("{} {}\n".format(p[0], p[1]))
    print(("Saved in {}".format(args.out_list)))


if __name__ == "__main__":
    main()
