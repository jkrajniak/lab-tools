import logging
import networkx as nx

from .md_libs import files_io

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def clean_atoms(atoms: list, monomers: list) -> dict:
    """Remove monomers from the list of atoms, make the atom ids continuous and return a dictionary with the mapping."""
    cleaned_atom_list = [a for a in atoms if a not in monomers]
    # make continous list
    mapping = {}
    for new_id, old_id in enumerate(cleaned_atom_list, 1):
        mapping[old_id] = new_id

    return mapping


def process(lmp: str, out: str):
    lammps_reader = files_io.LammpsReader()
    lammps_reader.read_data(lmp)
    logger.info(f"Read lammps data from {lmp}")
    logger.info(f"Number of atoms: {len(lammps_reader.atoms)}")

    lmp_graph = lammps_reader.get_simple_graph()
    logger.info(f"Created simple graph with {len(lmp_graph)} nodes and {len(lmp_graph.edges)} edges")

    # Find graph clusters
    clusters = list(nx.connected_components(lmp_graph))
    logger.info(f"Found {len(list(clusters))} clusters")
    # Get only the clusters with size 1
    monomers = [p for cluster in clusters for p in cluster if len(cluster) == 1]
    logger.info(f"Found {len(monomers)} monomers")

    # Remove monomers from the graph
    atom_mapping = clean_atoms(lammps_reader.atoms, monomers)

    # Remove atoms from lammps_reader
    lammps_reader.atoms = [atom_mapping[a] for a in lammps_reader.atoms if a not in monomers]
