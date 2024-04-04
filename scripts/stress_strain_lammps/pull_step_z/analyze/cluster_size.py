import argparse
import networkx as nx
import os
import numpy as np
import re
import cPickle

from md_libs import files_io

sys_ids = ("0_100", "25_75", "50_50", "75_25", "100_0")
total_number = 4000.0  # Total number of cross linked bonds

parser = argparse.ArgumentParser("Calculate cluster size")
parser.add_argument("sys_id", choices=sys_ids)

args = parser.parse_args()

# Prepare the mapping, atom_id -> residue_id
static_resid2mol_name = {x: "EPO" for x in range(1, 1801)}
static_resid2mol_name.update({x: "HDD" for x in range(1801, 2001)})
resid2mol_name = {}
for sys_id in sys_ids:
    resid2mol_name[sys_id] = static_resid2mol_name.copy()
    a, b = map(int, sys_id.split("_"))
    resid2mol_name[sys_id].update({x + 1: "JEF" for x in range(2000, 2000 + a * 10)})
    resid2mol_name[sys_id].update({x + 1: "IPD" for x in range(2000 + a * 10, 3000)})

a, b = map(int, args.sys_id.split("_"))
jef_range = (2000, 2001 + a * 10)
ipd_range = (2000 + a * 10, 3001)
print((jef_range, ipd_range))

output_data = []
re_filename = re.compile("data\.last\.[0-9.]+$")
cg_graphs = {}
for input_data in os.listdir("."):
    if re_filename.match(input_data):
        lr = files_io.LammpsReader(verbose=False)
        lr.read_data(input_data)

        g = nx.Graph()
        # g.add_nodes_from(range(1, 4001))
        for btype, bonds in lr.topology["bonds"].items():
            for b1, b2 in bonds:
                r1 = lr.atoms[b1]["res_id"]
                r2 = lr.atoms[b2]["res_id"]
                if r1 != r2:
                    g.add_edge(r1, r2, btype=btype)
        conversion = g.number_of_edges() / total_number
        connected_components = list(nx.connected_component_subgraphs(g))
        num_components = len(connected_components)
        component_size = [x.number_of_nodes() for x in connected_components]
        avg_deg_epo = np.average([d for k, d in g.degree().items() if k <= 1800])
        avg_deg_hdd = np.average([d for k, d in g.degree().items() if k > 1800 and k <= 2000])
        avg_deg_jef = np.average([d for k, d in g.degree().items() if k > jef_range[0] and k < jef_range[1]])
        avg_deg_ipd = np.average([d for k, d in g.degree().items() if k > ipd_range[0] and k < ipd_range[1]])
        if component_size:
            avg_component_size = np.average(component_size)
            largest_size = max(component_size)
        else:
            avg_component_size = 1.0
            largest_size = 1.0
        cg_graphs[conversion] = g
        print(
            (
                input_data,
                lr.timestep,
                conversion,
                num_components,
                largest_size,
                avg_component_size,
                avg_deg_epo,
                avg_deg_hdd,
                avg_deg_jef,
                avg_deg_ipd,
            )
        )
        output_data.append(
            [
                lr.timestep,
                conversion,
                num_components,
                largest_size,
                avg_component_size,
                avg_deg_epo,
                avg_deg_hdd,
                avg_deg_jef,
                avg_deg_ipd,
            ]
        )

output_data.sort(key=lambda x: x[0])
np.savetxt(
    "{}_cluster_data.csv".format(args.sys_id),
    np.nan_to_num(output_data),
    header="timestep conversion num_components largest_size avg_component_size avg_deg_epo avg_deg_hdd avg_deg_jef avg_deg_ipd",
)
with open("{}_graphs.pck".format(args.sys_id), "wb") as ob:
    cPickle.dump(cg_graphs, ob)
