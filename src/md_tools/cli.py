import click
import logging

from . import monomer_purger


@click.group()
@click.option('--debug/--no-debug', default=False)
def main(debug: bool):
    if debug:
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


@main.command(help='Purge monomers from LAMMPS network')
@click.option('-l', '--lmp', help='Input lammps configuration', required=True)
@click.option('-o', '--out', help='Output lammps configuration', required=True)
def purge(lmp: str, out: str):
    monomer_purger.process(lmp, out)


if __name__ == '__main__':
    main()
