import logging

from .md_libs import files_io

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def process(lmp: str, out: str):
    lammps_reader = files_io.LammpsReader()
    lammps_reader.read_data(lmp)
    logger.info(f"Read lammps data from {lmp}")
    logger.info(f"Number of atoms: {len(lammps_reader.atoms)}")

    lmp_graph = lammps_reader.get_graph()