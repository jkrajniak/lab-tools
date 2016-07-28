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
import re
from md_libs import files_io
import networkx as nx


def _args():
    parser = argparse.ArgumentParser('Generate exclusion lists')
    parser.add_argument('in_top', help='Input topology to process')
    parser.add_argument('out_list', help='Output exclusion list in the form of pairs')
    parser.add_argument('--nrexcl', type=int, default=3)
    parser.add_argument('--select_atoms', default=None, help=(
        'Generate exclusion list for subgroup of atoms. This is a comma separated'
        ' list of regular expressions (compiled by Python re module)'))

    return parser.parse_args()


def main():
    args = _args()
    input_top = files_io.GROMACSTopologyFile(args.in_top)
    input_top.read()

    print('Generate exclusion list with nrexcl={}'.format(args.nrexcl))

    selected_atom_ids = input_top.atoms.keys()
    if args.select_atoms:
        selected_atom_ids = []
        select_list = map(re.compile, args.select_atoms.split(','))
        for at_id, at_data in input_top.atoms.items():
            val = [re.match(m, at_data.name) for m in select_list]
            if any(val):
                selected_atom_ids.append(at_id)
    print('Selected {} atoms'.format(len(selected_atom_ids)))

    # Build a graph
    g = nx.Graph()
    for b1, b2 in input_top.bonds.keys():
        if b1 in selected_atom_ids and b2 in selected_atom_ids:
            g.add_edge(b1, b2)

    for b1, b2 in input_top.cross_bonds.keys():
        if b1 in selected_atom_ids and b2 in selected_atom_ids:
            g.add_edge(b1, b2)

    print('Generated graph, number of edges {}, nodes {}'.format(g.number_of_edges(), g.number_of_nodes()))
    all_paths = []
    connected_subgraphs = nx.connected_component_subgraphs(g)
    for sub_g in connected_subgraphs:
        node_list = sorted(sub_g.node)
        for i in range(sub_g.number_of_nodes()):
            for j in range(sub_g.number_of_nodes()):
                if i != j:
                    p = nx.all_simple_paths(sub_g, node_list[i], node_list[j], args.nrexcl)
                    all_paths.extend([tuple(sorted((x[0], x[-1]))) for x in p if x[0] != x[-1]])
    all_paths = sorted(set(all_paths))
    print('Generate {} exclusions'.format(len(all_paths)))

    with open(args.out_list, 'w') as output_file:
        for p in all_paths:
            output_file.write('{} {}\n'.format(p[0], p[1]))
    print('Saved in {}'.format(args.out_list))


if __name__ == '__main__':
    main()

