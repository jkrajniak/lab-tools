import dataclasses
import logging
import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)


def gen_bond_info(g: nx.Graph, natoms, ratio):
    group_a = int(natoms * ratio)
    group_b = natoms - group_a
    bond_info = []
    ntype = []
    for edge in g.edges:
        bond_info.append(edge)
        if edge[0] <= group_a and edge[1] <= group_a:
            ntype.append(1)
        elif edge[0] > group_a and edge[1] > group_a:
            ntype.append(2)
        else:
            ntype.append(3)
    return [ntype, bond_info]


def gen_bonded_tuples(g, num, bond_pair):
    """Generates tuples of different size, based on the graph and input edge.

    Args:
        g: The networkx Graph object.
        num: The length of the tuple.
        bond_pair: The edge which has to be included in all tuples.

    Returns:
        The set of all tuples of defined length from graph `g`.
    """
    b0, b1 = bond_pair
    paths = []
    if num > 3:
        for nb0 in g[b0]:
            paths.extend(nx.single_source_shortest_path(g, nb0, num - 1).values())
        for nb1 in g[b1]:
            paths.extend(nx.single_source_shortest_path(g, nb1, num - 1).values())

    paths.extend(nx.single_source_shortest_path(g, b0, num - 1).values())
    paths.extend(nx.single_source_shortest_path(g, b1, num - 1).values())
    output = set()
    for b in paths:
        if len(b) == num and b0 in b and b1 in b:
            if tuple(reversed(b)) not in output:
                output.add(tuple(b))
    return output


def gen_angle_info(g: nx.Graph, natoms: dict, ratio: float):
    group_a = int(natoms * ratio)
    group_b = natoms - group_a
    angle_info = []
    ntype = []
    for edge in g.edges:
        tem_angle = gen_bonded_tuples(g, 3, edge)
        for ijk in tem_angle:
            if tuple(reversed(ijk)) not in angle_info:
                angle_info.append(ijk)
    return [ntype, angle_info]


@dataclasses.dataclass
class LAMMPSWriterConfig:
    """The LAMMPS writer configuration object."""
    # The fraction of monomers in the system
    moo_frac: float
    # The scaling factor for the coordinates
    coo_scale: float
    # The length of the RMC box
    rmc_length: float


def write(self, lmp_cfg: LAMMPSWriterConfig, xzyfile_path: str, g: nx.Graph, out_path: str):
    """Generates the LAMMPS data file from the graph and the xzy file.

    Args:
        lmp_cfg: The LAMMPSWriterConfig object.
        xzyfile_path: The path to the xzy file.
        g: The networkx Graph object.
        out_path: The path to the output file.

    Returns:
        None
    """
    logger.info("Processing ATOM information ...")
    funcinfo = [g.degree(n) + 1 for n in g.nodes()]
    logger.info("Processing BOND information ...")
    natoms = len(g.nodes)
    bondinfo = gen_bond_info(g, natoms, lmp_cfg.mon_frac)
    print("Processing ANGLE information ...")
    angleinfo = gen_angle_info(g, natoms, lmp_cfg.mon_frac)
    totalAtoms = natoms
    groupA = int(natoms * lmp_cfg.mon_frac)
    coord = np.loadtxt(xzyfile_path)
    print("Writing LAMMPS data file ...")
    with open(out_path, 'w') as datainfo:
        datainfo.write('Polymer Network System\n')
        datainfo.write('\n')
        datainfo.write('{} atoms\n'.format(totalAtoms))
        datainfo.write('{} bonds\n'.format(len(bondinfo[0])))
        datainfo.write('{} angles\n'.format(len(angleinfo[0])))
        datainfo.write('\n')
        datainfo.write('5 atom types\n')
        datainfo.write('3 bond types\n')
        datainfo.write('4 angle types\n')
        datainfo.write('\n')
        datainfo.write('0.0000000 {} xlo xhi\n'.format(lmp_cfg.rmc_length * lmp_cfg.coo_scale))
        datainfo.write('0.0000000 {} ylo yhi\n'.format(lmp_cfg.rmc_length * lmp_cfg.coo_scale))
        datainfo.write('0.0000000 {} zlo zhi\n'.format(lmp_cfg.rmc_length * lmp_cfg.coo_scale))
        datainfo.write('\n')
        datainfo.write('Masses\n')
        datainfo.write('\n')
        datainfo.write('1 1.0\n')
        datainfo.write('2 1.0\n')
        datainfo.write('3 1.0\n')
        datainfo.write('4 1.0\n')
        datainfo.write('5 1.0\n')
        datainfo.write('\n')
        datainfo.write('Atoms\n')
        datainfo.write('\n')
        for i in range(totalAtoms):
            datainfo.write(
                '{:8d} {:8d} {:8d} {:10.6f} {:8.3f} {:8.3f} {:8.3f}\n'.format(i + 1, i + 1, funcinfo[i], 0.0,
                                                                              coord[:, 0][i] * lmp_cfg.coo_scale,
                                                                              coord[:, 1][i] * lmp_cfg.coo_scale,
                                                                              coord[:, 2][i] * lmp_cfg.coo_scale))
        datainfo.write('\n')
        datainfo.write('Bonds\n')
        datainfo.write('\n')
        for i in range(0, len(bondinfo[0])):
            datainfo.write(
                '{:8d} {:8d} {:8d} {:8d}\n'.format(i + 1, bondinfo[0][i], bondinfo[1][i][0], bondinfo[1][i][1]))
        datainfo.write('\n')
        datainfo.write('Angles\n')
        datainfo.write('\n')
        for i in range(0, len(angleinfo[0])):
            datainfo.write('{:8d} {:8d} {:8d} {:8d} {:8d}\n'.format(i + 1, angleinfo[0][i], angleinfo[1][i][0],
                                                                    angleinfo[1][i][1], angleinfo[1][i][2]))
        datainfo.write('\n')
    datainfo.close()
    print(len(bondinfo[0]), len(angleinfo[0]))
