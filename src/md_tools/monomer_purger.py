import logging
import networkx as nx

from .md_libs import files_io

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
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
    logger.info(f"Removing {len(monomers)} monomers from the lammps data")
    lammps_reader.atoms = {atom_mapping[a]: v for a, v in lammps_reader.atoms.items() if a not in monomers}
    lammps_reader.topology["bonds"] = {
        b: [(atom_mapping[b1], atom_mapping[b2]) for b1, b2 in blist]
            for b, blist in lammps_reader.topology["bonds"].items()
    }
    lammps_reader.topology["angles"] = {
        b: [(atom_mapping[b1], atom_mapping[b2], atom_mapping[b3]) for b1, b2, b3 in blist]
        for b, blist in lammps_reader.topology["angles"].items()
    }
    lammps_reader.topology["dihedrals"] = {
        b: [(atom_mapping[b1], atom_mapping[b2], atom_mapping[b3], atom_mapping[b4]) for b1, b2, b3, b4 in blist]
        for b, blist in lammps_reader.topology["dihedrals"].items()
    }
    lammps_reader.topology["impropers"] = {
        b: [(atom_mapping[b1], atom_mapping[b2], atom_mapping[b3], atom_mapping[b4]) for b1, b2, b3, b4 in blist]
        for b, blist in lammps_reader.topology["impropers"].items()
    }

    logger.info(f"Writing output to {out}")
    lammps_reader.write(out)
