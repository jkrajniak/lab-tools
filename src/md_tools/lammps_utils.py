import networkx as nx


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
            paths.extend(nx.single_source_shortest_path(g, nb0, num-1).values())
        for nb1 in g[b1]:
            paths.extend(nx.single_source_shortest_path(g, nb1, num-1).values())

    paths.extend(nx.single_source_shortest_path(g, b0, num-1).values())
    paths.extend(nx.single_source_shortest_path(g, b1, num-1).values())
    output = set()
    for b in paths:
        if len(b) == num and b0 in b and b1 in b:
            if tuple(reversed(b)) not in output:
                output.add(tuple(b))
    return output


def get_angle_info(g: nx.Graph, natoms: dict, ratio: float):
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