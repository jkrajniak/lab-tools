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
from matplotlib import pyplot as plt
import networkx as nx

from md_libs import files_io

parser = argparse.ArgumentParser(description="generator")
parser.add_argument("-t", "--top", help="GROMACS topology", required=True)
parser.add_argument("-d", "--draw", help="Draw graph", action="store_true", default=False)
parser.add_argument("-nocg", "--nocharges", help="Do not check charges", action="store_true", default=False)
parser.add_argument("-noang", "--noangles", help="Do not check angles", action="store_true", default=False)
parser.add_argument("-i", "--interactive", help="Run interactive shell", action="store_true", default=False)

args = parser.parse_args()

top = files_io.GROMACSTopologyFile(args.top)
top.read()

if not args.noangles or args.draw:
    g = nx.Graph()

    for at in list(top.atoms.values()):
        g.add_node(at.atom_id, name=at.name, chain_name=at.chain_name, chain_idx=at.chain_idx)

    for b in top.bonds:
        g.add_edge(*b)

color_node = {"O": "red", "C": "grey", "H": "white", "N": "lightblue"}


def get_color_nodes(gg, nodelist=None):
    if nodelist is None:
        nodelist = gg.nodes()
    return [color_node.get(gg.node[xx]["name"][0], "yellow") for xx in nodelist]


def draw_gg(gg, ommit=None):
    if ommit is None:
        ommit = []
    nodelist = [zz for zz, v in gg.node.items() if v["name"] not in ommit]
    nx.draw_graphviz(
        gg,
        with_labels=True,
        nodelist=nodelist,
        node_color=get_color_nodes(gg, nodelist),
        node_size=600,
        labels={xx: gg.node[xx]["name"] for xx in gg.nodes() if gg.node[xx]["name"] not in ommit},
    )
    plt.show()


if args.draw:
    color_nodes = [color_node.get(g.node[x]["name"][0], "yellow") for x in g.nodes()]

    nx.draw_graphviz(
        g,
        with_labels=True,
        node_color=color_nodes,
        node_size=600,
        labels={x: g.node[x]["name"] for x in g.nodes()},
    )
    plt.show()

if not args.noangles:
    # Generate angles
    angles = set([])
    dihedrals = set([])
    nodes = g.nodes()
    for idx, i in enumerate(nodes):
        for j in nodes[idx + 1 :]:
            path_i_j = nx.all_simple_paths(g, i, j, 4)
            for x in path_i_j:
                if len(x) == 3:
                    z, rev_z = tuple(x), tuple(reversed(x))
                    angles.add(z)
                    angles.add(rev_z)
                elif len(x) == 4:
                    z, rev_z = tuple(x), tuple(reversed(x))
                    dihedrals.add(z)
                    dihedrals.add(rev_z)

    # Compare with the topology itself.
    top_angles = set([a for a in top.angles])
    top_angles.update(set([tuple(reversed(a)) for a in top.angles]))
    print("Number of angles: top", len(top_angles) / 2, " and in the graph", len(angles) / 2)
    print("Set of angles equal", top_angles == angles)
    if not top_angles == angles:
        print("== Angles ==")
        print("Topology - graph", top_angles - angles)
        for idx, ang in enumerate((top_angles - angles)):
            print(idx, list(map(top.atoms.get, ang)))
        print("Graph - topology", angles - top_angles)
        for idx, ang in enumerate((angles - top_angles)):
            print(idx, list(map(top.atoms.get, ang)))

    print()
    # Check dihedrals
    top_dihedrals = set([d for d in top.dihedrals])
    top_dihedrals.update(set([tuple(reversed(d)) for d in top.dihedrals]))
    print("Number of dihedrals: top", len(top_dihedrals) / 2, " and in the graph", len(dihedrals) / 2)
    print("Set of dihedrals equal", top_dihedrals == dihedrals)
    if not top_dihedrals == dihedrals:
        print("== Dihedrals ==")
        print("Topology - graph", top_dihedrals - dihedrals)
        for idx, dih in enumerate((top_dihedrals - dihedrals)):
            print(idx, list(map(top.atoms.get, dih)))
        print("Graph - topology", dihedrals - top_dihedrals)
        for idx, dih in enumerate((dihedrals - top_dihedrals)):
            print(idx, [(x.atom_id, x.name, x.chain_name) for x in map(top.atoms.get, dih)])

if not args.nocharges:
    # Check if the charge is sum to 0
    total_charge = sum(x.charge if x.charge is not None else 0.0 for x in list(top.atoms.values()))
    print("Total charge:", total_charge)

if args.interactive:
    from IPython import embed

    embed()
