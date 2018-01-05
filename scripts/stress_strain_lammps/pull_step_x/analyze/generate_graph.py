
import sys
import networkx as nx

from md_libs import files_io

input_data = sys.argv[1]

lr = files_io.LammpsReader()
lr.read_data(input_data)

g = nx.Graph()

for at_id, at_data in lr.atoms.items():
    g.add_node(at_id, atom_type=at_data['atom_type'], res_id=at_data['res_id'])

for btype, bonds in lr.topology['bonds'].items():
    g.add_edges_from(bonds, btype=btype)

output_data = '{}_{:.2f}.pck'.format(input_data, (g.number_of_edges() - int(sys.argv[2]))/4000.0)
nx.write_gpickle(g, output_data)
print('Saved {}'.format(output_data))
