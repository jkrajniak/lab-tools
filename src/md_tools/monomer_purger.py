import logging
import networkx as nx

from .md_libs import files_io

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


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
    monomers = [cluster for cluster in clusters if len(cluster) == 1]
    logger.info(f"Found {len(monomers)} monomers")